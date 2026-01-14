#include "ExperimentController.h"

#include <QCoreApplication>
#include <QDateTime>
#include <QDir>
#include <QFileInfo>
#include <QRegularExpression>

namespace {
QString quoteArg(const QString& s) {
  if (s.isEmpty()) return "\"\"";
  if (!s.contains(' ') && !s.contains('\t') && !s.contains('"')) return s;
  QString escaped = s;
  escaped.replace('"', "\\\"");
  return "\"" + escaped + "\"";
}

QString joinCommand(const QString& program, const QStringList& args) {
  QStringList parts;
  parts.reserve(1 + args.size());
  parts << quoteArg(program);
  for (const auto& a : args) parts << quoteArg(a);
  return parts.join(' ');
}

double getDouble(const QVariantMap& m, const char* key, double fallback) {
  const auto it = m.find(QString::fromLatin1(key));
  return (it == m.end()) ? fallback : it->toDouble();
}

int getInt(const QVariantMap& m, const char* key, int fallback) {
  const auto it = m.find(QString::fromLatin1(key));
  return (it == m.end()) ? fallback : it->toInt();
}
}  // namespace

ExperimentController::ExperimentController(QObject* parent) : QObject(parent) {
  process_.setProcessChannelMode(QProcess::MergedChannels);
  connect(&process_, &QProcess::readyRead, this, &ExperimentController::handleProcessOutput);
  connect(
      &process_,
      qOverload<int, QProcess::ExitStatus>(&QProcess::finished),
      this,
      &ExperimentController::handleProcessFinished);
  connect(&process_, &QProcess::errorOccurred, this, &ExperimentController::handleProcessError);
}

RunQueueModel* ExperimentController::runModel() { return &model_; }

bool ExperimentController::running() const { return running_; }

int ExperimentController::totalRuns() const { return totalRuns_; }

int ExperimentController::finishedRuns() const { return finishedRuns_; }

int ExperimentController::selectedRun() const { return selectedRun_; }

void ExperimentController::setSelectedRun(int row) {
  if (selectedRun_ == row) return;
  selectedRun_ = row;
  emit selectedRunChanged();
  updateSelectedLog();
}

QString ExperimentController::selectedLog() const { return selectedLog_; }

QString ExperimentController::lastError() const { return lastError_; }

void ExperimentController::startRuns(const QVariantMap& config) {
  if (running_) {
    setLastError("已有任务在运行，请先停止。");
    return;
  }
  setLastError(QString());

  model_.clear();
  queue_.clear();
  totalRuns_ = 0;
  finishedRuns_ = 0;
  emit totalsChanged();

  const QVariantMap params = config.value("params").toMap();
  const QVariantMap sweeps = config.value("sweeps").toMap();

  QString simPath = normalizedPath(config.value("simPath").toString());
  if (simPath.isEmpty()) simPath = defaultSimPath();
  simPath = QFileInfo(simPath).absoluteFilePath();

  const QFileInfo simInfo(simPath);
  if (!simInfo.exists() || !simInfo.isFile()) {
    setLastError(QString("找不到实验程序：%1").arg(simPath));
    return;
  }

  const bool useMpi = config.value("useMpi").toBool();
  const int mpiRanks = std::max(1, config.value("mpiRanks").toInt());
  QString mpiLauncher = config.value("mpiLauncher").toString().trimmed();
  if (mpiLauncher.isEmpty()) mpiLauncher = "mpirun";

  QString outBase = normalizedPath(config.value("outBase").toString());
  if (outBase.isEmpty()) outBase = QDir::current().filePath("runs");
  const QString outBaseAbs = QDir(outBase).absolutePath();

  QString runPrefix = config.value("runPrefix").toString().trimmed();
  if (runPrefix.isEmpty()) runPrefix = "session";
  runPrefix = safeToken(runPrefix);

  const QString sessionStamp = QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss");
  const QString sessionDir = QDir(outBaseAbs).filePath(runPrefix + "_" + sessionStamp);
  if (!QDir().mkpath(sessionDir)) {
    setLastError(QString("无法创建输出目录：%1").arg(sessionDir));
    return;
  }

  const double u0 = getDouble(params, "u0", 0.05);
  const double con0 = getDouble(params, "con0", 0.2);
  const double sig = getDouble(params, "sig", 0.05);
  const double dt = getDouble(params, "dt", 0.05);
  const double dx = getDouble(params, "dx", 0.125);
  const int steps = std::max(1, getInt(params, "steps", 5002));
  const int mod = std::max(1, getInt(params, "mod", 25));
  const int grainx = std::max(1, getInt(params, "grainx", 80));
  const int grainy = std::max(1, getInt(params, "grainy", 80));
  const int grainz = std::max(1, getInt(params, "grainz", 1));
  const double axx = std::max(0.0, getDouble(params, "axx", 0.0));

  const QVariantMap u0Sweep = sweeps.value("u0").toMap();
  const QVector<double> u0Values = buildDoubleValues(
      u0Sweep.value("enabled").toBool(),
      u0,
      getDouble(u0Sweep, "start", u0),
      getDouble(u0Sweep, "end", u0),
      getDouble(u0Sweep, "step", 1.0));
  if (u0Values.isEmpty()) {
    setLastError("u0 扫描范围无效。");
    return;
  }

  const QVariantMap con0Sweep = sweeps.value("con0").toMap();
  const QVector<double> con0Values = buildDoubleValues(
      con0Sweep.value("enabled").toBool(),
      con0,
      getDouble(con0Sweep, "start", con0),
      getDouble(con0Sweep, "end", con0),
      getDouble(con0Sweep, "step", 1.0));
  if (con0Values.isEmpty()) {
    setLastError("con0 扫描范围无效。");
    return;
  }

  const QVariantMap stepsSweep = sweeps.value("steps").toMap();
  const QVector<int> stepsValues = buildIntValues(
      stepsSweep.value("enabled").toBool(),
      steps,
      getInt(stepsSweep, "start", steps),
      getInt(stepsSweep, "end", steps),
      getInt(stepsSweep, "step", 1));
  if (stepsValues.isEmpty()) {
    setLastError("steps 扫描范围无效。");
    return;
  }

  const qint64 runCount =
      static_cast<qint64>(u0Values.size()) * static_cast<qint64>(con0Values.size()) * static_cast<qint64>(stepsValues.size());
  if (runCount <= 0) {
    setLastError("组合数为 0。");
    return;
  }
  if (runCount > 5000) {
    setLastError(QString("组合数过大（%1），请加大步长或缩小范围。").arg(runCount));
    return;
  }

  int idx = 0;
  for (double u0v : u0Values) {
    for (double con0v : con0Values) {
      for (int stepsv : stepsValues) {
        ++idx;
        const QString name =
            QString("run_%1_u0_%2_con0_%3_steps_%4")
                .arg(idx, 4, 10, QChar('0'))
                .arg(safeToken(QString::number(u0v, 'g', 8)))
                .arg(safeToken(QString::number(con0v, 'g', 8)))
                .arg(stepsv);

        const QString outDir = QDir(sessionDir).filePath(name);
        QDir().mkpath(outDir);

        QStringList simArgs;
        simArgs << "--outdir" << outDir;
        simArgs << "--u0" << QString::number(u0v, 'g', 8);
        simArgs << "--con0" << QString::number(con0v, 'g', 8);
        simArgs << "--sig" << QString::number(sig, 'g', 8);
        simArgs << "--dt" << QString::number(dt, 'g', 8);
        simArgs << "--dx" << QString::number(dx, 'g', 8);
        simArgs << "--steps" << QString::number(stepsv);
        simArgs << "--mod" << QString::number(mod);
        simArgs << "--grainx" << QString::number(grainx);
        simArgs << "--grainy" << QString::number(grainy);
        simArgs << "--grainz" << QString::number(grainz);
        simArgs << "--axx" << QString::number(axx, 'g', 8);

        RunSpec spec;
        spec.name = name;
        spec.outDir = outDir;
        spec.guiLogPath = QDir(outDir).filePath("gui_process.log");

        if (useMpi) {
          spec.program = mpiLauncher;
          spec.args << "-np" << QString::number(mpiRanks);
          spec.args << simPath;
          spec.args << simArgs;
        } else {
          spec.program = simPath;
          spec.args = simArgs;
        }

        spec.commandString = joinCommand(spec.program, spec.args);

        RunQueueModel::Item item;
        item.name = name;
        item.status = "Pending";
        item.outDir = outDir;
        item.command = spec.commandString;
        spec.row = model_.addItem(item);

        queue_.enqueue(spec);
      }
    }
  }

  totalRuns_ = idx;
  finishedRuns_ = 0;
  emit totalsChanged();

  setRunning(true);
  startNext();
}

void ExperimentController::stop() {
  if (!running_) return;
  queue_.clear();
  if (process_.state() != QProcess::NotRunning) {
    process_.kill();
    process_.waitForFinished(2000);
  }
  if (currentRow_ >= 0) {
    model_.setStatus(currentRow_, "Stopped");
    model_.markFinished(currentRow_, -1);
    ++finishedRuns_;
    emit totalsChanged();
  }
  currentRow_ = -1;
  if (currentGuiLog_.isOpen()) currentGuiLog_.close();
  setRunning(false);
}

void ExperimentController::clear() {
  stop();
  model_.clear();
  selectedRun_ = -1;
  selectedLog_.clear();
  emit selectedRunChanged();
  emit selectedLogChanged();
}

void ExperimentController::setRunning(bool v) {
  if (running_ == v) return;
  running_ = v;
  emit runningChanged();
}

void ExperimentController::setLastError(const QString& err) {
  if (lastError_ == err) return;
  lastError_ = err;
  emit lastErrorChanged();
}

void ExperimentController::updateSelectedLog() {
  const auto* item = model_.itemAt(selectedRun_);
  const QString next = item ? item->logTail : QString();
  if (selectedLog_ == next) return;
  selectedLog_ = next;
  emit selectedLogChanged();
}

void ExperimentController::startNext() {
  if (process_.state() != QProcess::NotRunning) return;
  if (queue_.isEmpty()) {
    currentRow_ = -1;
    setRunning(false);
    return;
  }

  const RunSpec spec = queue_.dequeue();
  currentRow_ = spec.row;

  model_.setStatus(currentRow_, "Running");
  model_.markStarted(currentRow_);

  if (currentGuiLog_.isOpen()) currentGuiLog_.close();
  currentGuiLog_.setFileName(spec.guiLogPath);
  if (!currentGuiLog_.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Text)) {
    model_.appendLog(currentRow_, QString("Warning: Failed to open %1\n").arg(spec.guiLogPath));
  }

  process_.setProgram(spec.program);
  process_.setArguments(spec.args);
  process_.setWorkingDirectory(spec.outDir);
  process_.start();
}

void ExperimentController::handleProcessOutput() {
  if (currentRow_ < 0) return;
  const QByteArray bytes = process_.readAll();
  if (bytes.isEmpty()) return;

  const QString text = QString::fromLocal8Bit(bytes);
  model_.appendLog(currentRow_, text);
  if (currentGuiLog_.isOpen()) {
    currentGuiLog_.write(bytes);
    currentGuiLog_.flush();
  }
  if (selectedRun_ == currentRow_) updateSelectedLog();
}

void ExperimentController::handleProcessFinished(int exitCode, QProcess::ExitStatus exitStatus) {
  if (currentGuiLog_.isOpen()) currentGuiLog_.close();
  if (currentRow_ >= 0) {
    const bool ok = (exitStatus == QProcess::NormalExit && exitCode == 0);
    model_.setStatus(currentRow_, ok ? "Done" : "Failed");
    model_.markFinished(currentRow_, exitCode);
    ++finishedRuns_;
    emit totalsChanged();
    if (selectedRun_ == currentRow_) updateSelectedLog();
  }
  currentRow_ = -1;
  startNext();
}

void ExperimentController::handleProcessError(QProcess::ProcessError error) {
  Q_UNUSED(error);
  if (currentRow_ >= 0) {
    model_.setStatus(currentRow_, "Error");
    model_.markFinished(currentRow_, -1);
    ++finishedRuns_;
    emit totalsChanged();
    if (selectedRun_ == currentRow_) updateSelectedLog();
  }
  currentRow_ = -1;
  startNext();
}

QString ExperimentController::defaultSimPath() {
  const QString dir = QCoreApplication::applicationDirPath();
#ifdef Q_OS_WIN
  return QDir(dir).filePath("misfit_sim.exe");
#else
  return QDir(dir).filePath("misfit_sim");
#endif
}

QString ExperimentController::normalizedPath(QString p) {
  return QDir::cleanPath(p.trimmed());
}

QString ExperimentController::safeToken(QString s) {
  s = s.trimmed();
  s.replace(QRegularExpression("[^A-Za-z0-9._-]+"), "_");
  s.replace(QRegularExpression("_+"), "_");
  s.replace(QRegularExpression("^_+|_+$"), "");
  if (s.isEmpty()) s = "v";
  return s;
}

QVector<double> ExperimentController::buildDoubleValues(bool enabled, double single, double start, double end, double step) {
  if (!enabled) return {single};
  if (!(step > 0.0)) return {};
  if (end < start) return {};
  QVector<double> values;
  for (int i = 0;; ++i) {
    const double v = start + static_cast<double>(i) * step;
    if (v > end + 1e-12) break;
    values.push_back(v);
    if (values.size() > 5000) break;
  }
  return values;
}

QVector<int> ExperimentController::buildIntValues(bool enabled, int single, int start, int end, int step) {
  if (!enabled) return {single};
  if (step <= 0) return {};
  if (end < start) return {};
  QVector<int> values;
  for (int v = start; v <= end; v += step) {
    values.push_back(v);
    if (values.size() > 5000) break;
  }
  return values;
}

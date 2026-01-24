#include "MainWindow.h"

#include <QCheckBox>
#include <QCoreApplication>
#include <QDateTime>
#include <QDesktopServices>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QSpinBox>
#include <QStandardPaths>
#include <QStyle>
#include <QTabWidget>
#include <QToolButton>
#include <QUrl>
#include <QVBoxLayout>
#include <QFontDatabase>

static QDoubleSpinBox* makeDoubleSpin(double min, double max, double value, int decimals, double step) {
    auto* box = new QDoubleSpinBox;
    box->setRange(min, max);
    box->setDecimals(decimals);
    box->setSingleStep(step);
    box->setValue(value);
    box->setKeyboardTracking(false);
    return box;
}

static QSpinBox* makeIntSpin(int min, int max, int value, int step) {
    auto* box = new QSpinBox;
    box->setRange(min, max);
    box->setSingleStep(step);
    box->setValue(value);
    box->setKeyboardTracking(false);
    return box;
}

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("PFC-Exp-GUI");
    resize(1100, 780);

    auto* tabs = new QTabWidget;
    auto* singleTab = new QWidget;
    auto* batchTab = new QWidget;
    tabs->addTab(singleTab, "Single Run");
    tabs->addTab(batchTab, "Batch Sweep");

    log_ = new QPlainTextEdit;
    log_->setReadOnly(true);
    log_->setMaximumBlockCount(5000);
    {
        QFont f = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        if (f.pointSizeF() > 0) {
            f.setPointSizeF(std::max(9.0, f.pointSizeF()));
        } else {
            f.setPointSize(10);
        }
        log_->setFont(f);
    }

    auto* logBox = new QGroupBox("Log");
    auto* logLayout = new QVBoxLayout;
    logLayout->addWidget(log_);
    logLayout->setContentsMargins(0, 0, 0, 0);
    logBox->setLayout(logLayout);

    auto* central = new QWidget;
    auto* centralLayout = new QVBoxLayout;
    centralLayout->setContentsMargins(14, 14, 14, 14);
    centralLayout->setSpacing(12);
    centralLayout->addWidget(tabs);
    centralLayout->addWidget(logBox, 1);
    central->setLayout(centralLayout);
    setCentralWidget(central);

    // ---------- Single tab ----------
    auto* singleLayout = new QVBoxLayout;
    singleLayout->setSpacing(12);
    singleTab->setLayout(singleLayout);

    auto* singleParamsBox = new QGroupBox("Parameters");
    auto* singleParamsForm = new QFormLayout;
    singleParamsForm->setVerticalSpacing(8);
    singleParamsForm->setHorizontalSpacing(12);

    singleU0_ = makeDoubleSpin(-2.0, 2.0, 0.05, 6, 0.01);
    singleCon0_ = makeDoubleSpin(0.0, 1.0, 0.2, 6, 0.01);
    singleSig_ = makeDoubleSpin(0.0, 2.0, 0.05, 6, 0.01);
    singleDt_ = makeDoubleSpin(1e-6, 1.0, 0.05, 6, 0.01);
    singleDx_ = makeDoubleSpin(1e-6, 10.0, 0.125, 6, 0.01);
    singleSteps_ = makeIntSpin(1, 100000000, 200, 10);
    singleMod_ = makeIntSpin(1, 100000000, 25, 1);
    singleSeed_ = makeIntSpin(0, 2147483647, 20200604, 1);
    singleGrainX_ = makeIntSpin(1, 4096, 16, 1);
    singleGrainY_ = makeIntSpin(1, 4096, 16, 1);
    singleGrainZ_ = makeIntSpin(1, 4096, 1, 1);
    singleAxx_ = makeDoubleSpin(0.0, 1000.0, 0.0, 6, 0.1);

    singleParamsForm->addRow("u0 (density mean)", singleU0_);
    singleParamsForm->addRow("con0 (concentration mean)", singleCon0_);
    singleParamsForm->addRow("sig", singleSig_);
    singleParamsForm->addRow("dt", singleDt_);
    singleParamsForm->addRow("dx", singleDx_);
    singleParamsForm->addRow("steps", singleSteps_);
    singleParamsForm->addRow("mod (checkpoint interval)", singleMod_);
    singleParamsForm->addRow("seed (rng)", singleSeed_);
    singleParamsForm->addRow("grainx", singleGrainX_);
    singleParamsForm->addRow("grainy", singleGrainY_);
    singleParamsForm->addRow("grainz", singleGrainZ_);
    singleParamsForm->addRow("axx (noise intensity)", singleAxx_);
    singleParamsBox->setLayout(singleParamsForm);

    auto* singleOutBox = new QGroupBox("Output");
    auto* singleOutForm = new QFormLayout;
    singleBaseOutDir_ = new QLineEdit(defaultBaseOutputDir());
    singleRunName_ = new QLineEdit;
    singleRunName_->setPlaceholderText("Optional, e.g. test_run_01 (default: timestamp)");

    auto* singleBrowse = new QToolButton;
    singleBrowse->setText("Browse…");
    singleBrowse->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    singleBrowse->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    connect(singleBrowse, &QToolButton::clicked, this, &MainWindow::browseSingleOutDir);

    auto* baseOutRow = new QWidget;
    auto* baseOutRowLayout = new QHBoxLayout;
    baseOutRowLayout->setContentsMargins(0, 0, 0, 0);
    baseOutRowLayout->addWidget(singleBaseOutDir_, 1);
    baseOutRowLayout->addWidget(singleBrowse);
    baseOutRow->setLayout(baseOutRowLayout);

    singleOutForm->addRow("Base output dir", baseOutRow);
    singleOutForm->addRow("Run folder name", singleRunName_);
    singleOutBox->setLayout(singleOutForm);

    auto* singleButtons = new QWidget;
    auto* singleButtonsLayout = new QHBoxLayout;
    singleButtonsLayout->setContentsMargins(0, 0, 0, 0);
    auto* singleRun = new QPushButton("Run");
    auto* singleStop = new QPushButton("Stop");
    auto* openOut = new QPushButton("Open Output");
    singleRun->setProperty("primary", true);
    singleStop->setProperty("danger", true);
    singleRun->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    singleStop->setIcon(style()->standardIcon(QStyle::SP_MediaStop));
    openOut->setIcon(style()->standardIcon(QStyle::SP_DirOpenIcon));
    connect(singleRun, &QPushButton::clicked, this, &MainWindow::startSingleRun);
    connect(singleStop, &QPushButton::clicked, this, &MainWindow::stopRun);
    connect(openOut, &QPushButton::clicked, this, &MainWindow::openCurrentOutputDir);
    singleButtonsLayout->addWidget(singleRun);
    singleButtonsLayout->addWidget(singleStop);
    singleButtonsLayout->addStretch(1);
    singleButtonsLayout->addWidget(openOut);
    singleButtons->setLayout(singleButtonsLayout);

    singleLayout->addWidget(singleParamsBox);
    singleLayout->addWidget(singleOutBox);
    singleLayout->addWidget(singleButtons);
    singleLayout->addStretch(1);

    // ---------- Batch tab ----------
    auto* batchLayout = new QVBoxLayout;
    batchLayout->setSpacing(12);
    batchTab->setLayout(batchLayout);

    auto* batchFixedBox = new QGroupBox("Fixed Parameters");
    auto* batchFixedForm = new QFormLayout;
    batchFixedForm->setVerticalSpacing(8);
    batchFixedForm->setHorizontalSpacing(12);
    batchSig_ = makeDoubleSpin(0.0, 2.0, 0.05, 6, 0.01);
    batchDt_ = makeDoubleSpin(1e-6, 1.0, 0.05, 6, 0.01);
    batchDx_ = makeDoubleSpin(1e-6, 10.0, 0.125, 6, 0.01);
    batchMod_ = makeIntSpin(1, 100000000, 25, 1);
    batchSeed_ = makeIntSpin(0, 2147483647, 20200604, 1);
    batchGrainX_ = makeIntSpin(1, 4096, 16, 1);
    batchGrainY_ = makeIntSpin(1, 4096, 16, 1);
    batchGrainZ_ = makeIntSpin(1, 4096, 1, 1);
    batchAxx_ = makeDoubleSpin(0.0, 1000.0, 0.0, 6, 0.1);

    batchFixedForm->addRow("sig", batchSig_);
    batchFixedForm->addRow("dt", batchDt_);
    batchFixedForm->addRow("dx", batchDx_);
    batchFixedForm->addRow("mod (checkpoint interval)", batchMod_);
    batchFixedForm->addRow("seed (rng)", batchSeed_);
    batchFixedForm->addRow("grainx", batchGrainX_);
    batchFixedForm->addRow("grainy", batchGrainY_);
    batchFixedForm->addRow("grainz", batchGrainZ_);
    batchFixedForm->addRow("axx (noise intensity)", batchAxx_);
    batchFixedBox->setLayout(batchFixedForm);

    auto* batchSweepBox = new QGroupBox("Sweep Parameters (Cartesian product)");
    auto* sweepForm = new QFormLayout;
    sweepForm->setVerticalSpacing(8);
    sweepForm->setHorizontalSpacing(12);

    auto* u0Row = new QWidget;
    auto* u0RowLayout = new QHBoxLayout;
    u0RowLayout->setContentsMargins(0, 0, 0, 0);
    sweepU0Enable_ = new QCheckBox("Sweep");
    sweepU0Single_ = makeDoubleSpin(-2.0, 2.0, 0.05, 6, 0.01);
    sweepU0Start_ = makeDoubleSpin(-2.0, 2.0, 0.0, 6, 0.01);
    sweepU0End_ = makeDoubleSpin(-2.0, 2.0, 0.2, 6, 0.01);
    sweepU0Step_ = makeDoubleSpin(1e-9, 10.0, 0.05, 6, 0.01);
    u0RowLayout->addWidget(sweepU0Enable_);
    u0RowLayout->addWidget(new QLabel("Single"));
    u0RowLayout->addWidget(sweepU0Single_);
    u0RowLayout->addSpacing(10);
    u0RowLayout->addWidget(new QLabel("Start"));
    u0RowLayout->addWidget(sweepU0Start_);
    u0RowLayout->addWidget(new QLabel("End"));
    u0RowLayout->addWidget(sweepU0End_);
    u0RowLayout->addWidget(new QLabel("Step"));
    u0RowLayout->addWidget(sweepU0Step_);
    u0Row->setLayout(u0RowLayout);
    sweepForm->addRow("u0", u0Row);

    auto* con0Row = new QWidget;
    auto* con0RowLayout = new QHBoxLayout;
    con0RowLayout->setContentsMargins(0, 0, 0, 0);
    sweepCon0Enable_ = new QCheckBox("Sweep");
    sweepCon0Single_ = makeDoubleSpin(0.0, 1.0, 0.2, 6, 0.01);
    sweepCon0Start_ = makeDoubleSpin(0.0, 1.0, 0.1, 6, 0.01);
    sweepCon0End_ = makeDoubleSpin(0.0, 1.0, 0.3, 6, 0.01);
    sweepCon0Step_ = makeDoubleSpin(1e-9, 1.0, 0.05, 6, 0.01);
    con0RowLayout->addWidget(sweepCon0Enable_);
    con0RowLayout->addWidget(new QLabel("Single"));
    con0RowLayout->addWidget(sweepCon0Single_);
    con0RowLayout->addSpacing(10);
    con0RowLayout->addWidget(new QLabel("Start"));
    con0RowLayout->addWidget(sweepCon0Start_);
    con0RowLayout->addWidget(new QLabel("End"));
    con0RowLayout->addWidget(sweepCon0End_);
    con0RowLayout->addWidget(new QLabel("Step"));
    con0RowLayout->addWidget(sweepCon0Step_);
    con0Row->setLayout(con0RowLayout);
    sweepForm->addRow("con0", con0Row);

    auto* stepsRow = new QWidget;
    auto* stepsRowLayout = new QHBoxLayout;
    stepsRowLayout->setContentsMargins(0, 0, 0, 0);
    sweepStepsEnable_ = new QCheckBox("Sweep");
    sweepStepsSingle_ = makeIntSpin(1, 100000000, 200, 10);
    sweepStepsStart_ = makeIntSpin(1, 100000000, 100, 10);
    sweepStepsEnd_ = makeIntSpin(1, 100000000, 300, 10);
    sweepStepsStep_ = makeIntSpin(1, 100000000, 100, 10);
    stepsRowLayout->addWidget(sweepStepsEnable_);
    stepsRowLayout->addWidget(new QLabel("Single"));
    stepsRowLayout->addWidget(sweepStepsSingle_);
    stepsRowLayout->addSpacing(10);
    stepsRowLayout->addWidget(new QLabel("Start"));
    stepsRowLayout->addWidget(sweepStepsStart_);
    stepsRowLayout->addWidget(new QLabel("End"));
    stepsRowLayout->addWidget(sweepStepsEnd_);
    stepsRowLayout->addWidget(new QLabel("Step"));
    stepsRowLayout->addWidget(sweepStepsStep_);
    stepsRow->setLayout(stepsRowLayout);
    sweepForm->addRow("steps", stepsRow);

    batchSweepBox->setLayout(sweepForm);

    auto* batchOutBox = new QGroupBox("Output");
    auto* batchOutForm = new QFormLayout;
    batchBaseOutDir_ = new QLineEdit(defaultBaseOutputDir());
    batchName_ = new QLineEdit;
    batchName_->setPlaceholderText("Optional (default: timestamp)");
    auto* batchBrowse = new QToolButton;
    batchBrowse->setText("Browse…");
    batchBrowse->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    batchBrowse->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    connect(batchBrowse, &QToolButton::clicked, this, &MainWindow::browseBatchOutDir);

    auto* batchBaseRow = new QWidget;
    auto* batchBaseRowLayout = new QHBoxLayout;
    batchBaseRowLayout->setContentsMargins(0, 0, 0, 0);
    batchBaseRowLayout->addWidget(batchBaseOutDir_, 1);
    batchBaseRowLayout->addWidget(batchBrowse);
    batchBaseRow->setLayout(batchBaseRowLayout);

    batchOutForm->addRow("Base output dir", batchBaseRow);
    batchOutForm->addRow("Batch folder name", batchName_);
    batchOutBox->setLayout(batchOutForm);

    batchPreview_ = new QLabel;
    batchPreview_->setProperty("hint", true);
    batchPreview_->setWordWrap(true);
    batchProgress_ = new QProgressBar;
    batchProgress_->setRange(0, 100);
    batchProgress_->setValue(0);

    auto* batchButtons = new QWidget;
    auto* batchButtonsLayout = new QHBoxLayout;
    batchButtonsLayout->setContentsMargins(0, 0, 0, 0);
    auto* batchRun = new QPushButton("Run Batch");
    auto* batchStop = new QPushButton("Stop");
    auto* openOut2 = new QPushButton("Open Output");
    batchRun->setProperty("primary", true);
    batchStop->setProperty("danger", true);
    batchRun->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    batchStop->setIcon(style()->standardIcon(QStyle::SP_MediaStop));
    openOut2->setIcon(style()->standardIcon(QStyle::SP_DirOpenIcon));
    connect(batchRun, &QPushButton::clicked, this, &MainWindow::startBatchRun);
    connect(batchStop, &QPushButton::clicked, this, &MainWindow::stopRun);
    connect(openOut2, &QPushButton::clicked, this, &MainWindow::openCurrentOutputDir);
    batchButtonsLayout->addWidget(batchRun);
    batchButtonsLayout->addWidget(batchStop);
    batchButtonsLayout->addStretch(1);
    batchButtonsLayout->addWidget(openOut2);
    batchButtons->setLayout(batchButtonsLayout);

    batchLayout->addWidget(batchFixedBox);
    batchLayout->addWidget(batchSweepBox);
    batchLayout->addWidget(batchOutBox);
    batchLayout->addWidget(batchPreview_);
    batchLayout->addWidget(batchProgress_);
    batchLayout->addWidget(batchButtons);
    batchLayout->addStretch(1);

    auto connectValueChanged = [&](QObject* widget) {
        if (auto* d = qobject_cast<QDoubleSpinBox*>(widget)) {
            connect(d, &QDoubleSpinBox::valueChanged, this, &MainWindow::updateBatchPreview);
        } else if (auto* s = qobject_cast<QSpinBox*>(widget)) {
            connect(s, &QSpinBox::valueChanged, this, &MainWindow::updateBatchPreview);
        } else if (auto* c = qobject_cast<QCheckBox*>(widget)) {
            connect(c, &QCheckBox::toggled, this, &MainWindow::updateBatchPreview);
        } else if (auto* e = qobject_cast<QLineEdit*>(widget)) {
            connect(e, &QLineEdit::textChanged, this, &MainWindow::updateBatchPreview);
        }
    };

    connectValueChanged(sweepU0Enable_);
    connectValueChanged(sweepU0Single_);
    connectValueChanged(sweepU0Start_);
    connectValueChanged(sweepU0End_);
    connectValueChanged(sweepU0Step_);
    connectValueChanged(sweepCon0Enable_);
    connectValueChanged(sweepCon0Single_);
    connectValueChanged(sweepCon0Start_);
    connectValueChanged(sweepCon0End_);
    connectValueChanged(sweepCon0Step_);
    connectValueChanged(sweepStepsEnable_);
    connectValueChanged(sweepStepsSingle_);
    connectValueChanged(sweepStepsStart_);
    connectValueChanged(sweepStepsEnd_);
    connectValueChanged(sweepStepsStep_);
    connectValueChanged(batchName_);
    connectValueChanged(batchBaseOutDir_);

    updateBatchPreview();
    setRunningUi(false);
}

QString MainWindow::cliPath() const {
#ifdef Q_OS_WIN
    const QString cliName = "pfc-exp-cli.exe";
#else
    const QString cliName = "pfc-exp-cli";
#endif
    const QString dir = QCoreApplication::applicationDirPath();
    const QString candidate = QDir(dir).filePath(cliName);
    if (QFileInfo::exists(candidate)) return candidate;
    return candidate;
}

QString MainWindow::defaultBaseOutputDir() const {
    const QString docs = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    if (!docs.isEmpty()) {
        return QDir(docs).filePath("PFC-Exp-GUI/outputs");
    }
    return QDir::current().filePath("outputs");
}

bool MainWindow::ensureDir(const QString& dirPath, QString* errorOut) const {
    if (dirPath.trimmed().isEmpty()) {
        if (errorOut) *errorOut = "Output directory is empty.";
        return false;
    }
    QDir dir(dirPath);
    if (dir.exists()) return true;
    if (!dir.mkpath(".")) {
        if (errorOut) *errorOut = "Failed to create output directory: " + dirPath;
        return false;
    }
    return true;
}

QString MainWindow::makeTimestampedDirName(const QString& prefix) const {
    const QString ts = QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss");
    return prefix + "_" + ts;
}

QString MainWindow::formatRunLabel(int index, int total) const {
    if (total <= 1) return "Run";
    return QString("Run %1 / %2").arg(index).arg(total);
}

RunParams MainWindow::singleParams() const {
    RunParams p;
    p.u0 = singleU0_->value();
    p.con0 = singleCon0_->value();
    p.sig = singleSig_->value();
    p.dt = singleDt_->value();
    p.dx = singleDx_->value();
    p.steps = singleSteps_->value();
    p.mod = singleMod_->value();
    p.seed = singleSeed_->value();
    p.grainx = singleGrainX_->value();
    p.grainy = singleGrainY_->value();
    p.grainz = singleGrainZ_->value();
    p.axx = singleAxx_->value();
    return p;
}

RunParams MainWindow::batchFixedParams() const {
    RunParams p;
    p.sig = batchSig_->value();
    p.dt = batchDt_->value();
    p.dx = batchDx_->value();
    p.mod = batchMod_->value();
    p.seed = batchSeed_->value();
    p.grainx = batchGrainX_->value();
    p.grainy = batchGrainY_->value();
    p.grainz = batchGrainZ_->value();
    p.axx = batchAxx_->value();
    return p;
}

QStringList MainWindow::buildArgs(const RunParams& params, const QString& outDir) const {
    auto f = [](double v) { return QString::number(v, 'g', 10); };
    QStringList args;
    args << "--u0" << f(params.u0);
    args << "--con0" << f(params.con0);
    args << "--sig" << f(params.sig);
    args << "--dt" << f(params.dt);
    args << "--dx" << f(params.dx);
    args << "--steps" << QString::number(params.steps);
    args << "--mod" << QString::number(params.mod);
    args << "--seed" << QString::number(params.seed);
    args << "--grainx" << QString::number(params.grainx);
    args << "--grainy" << QString::number(params.grainy);
    args << "--grainz" << QString::number(params.grainz);
    args << "--axx" << f(params.axx);
    args << "--outdir" << outDir;
    return args;
}

bool MainWindow::writeParamsJson(const RunJob& job, QString* errorOut) const {
    QJsonObject obj;
    obj["mode"] = job.mode;
    obj["run_index"] = job.index;
    obj["run_total"] = job.total;
    obj["out_dir"] = job.outDir;

    QJsonObject p;
    p["u0"] = job.params.u0;
    p["con0"] = job.params.con0;
    p["sig"] = job.params.sig;
    p["dt"] = job.params.dt;
    p["dx"] = job.params.dx;
    p["steps"] = job.params.steps;
    p["mod"] = job.params.mod;
    p["seed"] = job.params.seed;
    p["grainx"] = job.params.grainx;
    p["grainy"] = job.params.grainy;
    p["grainz"] = job.params.grainz;
    p["axx"] = job.params.axx;
    obj["params"] = p;

    const QString path = QDir(job.outDir).filePath("params.json");
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (errorOut) *errorOut = "Failed to write params.json: " + path;
        return false;
    }
    file.write(QJsonDocument(obj).toJson(QJsonDocument::Indented));
    file.close();
    return true;
}

static QVector<double> makeDoubleRange(double start, double end, double step, QString* err) {
    QVector<double> out;
    if (step <= 0.0) {
        if (err) *err = "Step must be > 0.";
        return out;
    }
    if (end < start) {
        if (err) *err = "End must be >= start.";
        return out;
    }
    const int maxCount = 10000;
    for (int i = 0; i < maxCount; ++i) {
        const double v = start + step * i;
        if (v > end + 1e-12) break;
        out.push_back(v);
    }
    if (out.isEmpty()) out.push_back(start);
    return out;
}

static QVector<int> makeIntRange(int start, int end, int step, QString* err) {
    QVector<int> out;
    if (step <= 0) {
        if (err) *err = "Step must be > 0.";
        return out;
    }
    if (end < start) {
        if (err) *err = "End must be >= start.";
        return out;
    }
    const int maxCount = 100000;
    for (int v = start, i = 0; v <= end && i < maxCount; v += step, ++i) out.push_back(v);
    if (out.isEmpty()) out.push_back(start);
    return out;
}

QVector<RunJob> MainWindow::buildBatchJobs(const QString& baseDir, const QString& batchName, QString* errorOut) const {
    const int maxRuns = 5000;

    const RunParams fixed = batchFixedParams();

    QString err;
    const QVector<double> u0Values = sweepU0Enable_->isChecked()
        ? makeDoubleRange(sweepU0Start_->value(), sweepU0End_->value(), sweepU0Step_->value(), &err)
        : QVector<double>{sweepU0Single_->value()};
    if (!err.isEmpty()) {
        if (errorOut) *errorOut = "u0: " + err;
        return {};
    }

    const QVector<double> con0Values = sweepCon0Enable_->isChecked()
        ? makeDoubleRange(sweepCon0Start_->value(), sweepCon0End_->value(), sweepCon0Step_->value(), &err)
        : QVector<double>{sweepCon0Single_->value()};
    if (!err.isEmpty()) {
        if (errorOut) *errorOut = "con0: " + err;
        return {};
    }

    const QVector<int> stepsValues = sweepStepsEnable_->isChecked()
        ? makeIntRange(sweepStepsStart_->value(), sweepStepsEnd_->value(), sweepStepsStep_->value(), &err)
        : QVector<int>{sweepStepsSingle_->value()};
    if (!err.isEmpty()) {
        if (errorOut) *errorOut = "steps: " + err;
        return {};
    }

    const long long total =
        static_cast<long long>(u0Values.size()) * static_cast<long long>(con0Values.size()) * static_cast<long long>(stepsValues.size());
    if (total <= 0) {
        if (errorOut) *errorOut = "No runs to execute.";
        return {};
    }
    if (total > maxRuns) {
        if (errorOut) *errorOut = QString("Too many runs (%1). Reduce ranges or steps (limit: %2).").arg(total).arg(maxRuns);
        return {};
    }
    if (baseDir.trimmed().isEmpty()) {
        if (errorOut) *errorOut = "Base output dir is empty.";
        return {};
    }
    if (batchName.trimmed().isEmpty()) {
        if (errorOut) *errorOut = "Batch folder name is empty.";
        return {};
    }
    const QString batchRoot = QDir(baseDir).filePath(batchName);

    QVector<RunJob> jobs;
    jobs.reserve(static_cast<int>(total));
    int idx = 1;
    for (const double u0 : u0Values) {
        for (const double con0 : con0Values) {
            for (const int steps : stepsValues) {
                RunJob job;
                job.mode = "batch";
                job.index = idx;
                job.total = static_cast<int>(total);
                job.params = fixed;
                job.params.u0 = u0;
                job.params.con0 = con0;
                job.params.steps = steps;
                job.outDir = QDir(batchRoot).filePath(QString("run_%1").arg(idx, 4, 10, QChar('0')));
                jobs.push_back(job);
                ++idx;
            }
        }
    }
    return jobs;
}

void MainWindow::appendLog(const QString& text) {
    log_->appendPlainText(text);
}

void MainWindow::setRunningUi(bool running) {
    (void)running;
}

void MainWindow::launchJob(const RunJob& job) {
    QString err;
    if (!ensureDir(job.outDir, &err)) {
        QMessageBox::critical(this, "Output error", err);
        return;
    }
    if (!writeParamsJson(job, &err)) {
        QMessageBox::critical(this, "Output error", err);
        return;
    }

    const QString program = cliPath();
    if (!QFileInfo::exists(program)) {
        QMessageBox::critical(this, "CLI not found",
                              "Cannot find experiment CLI next to the GUI executable:\n" + program +
                                  "\n\nMake sure `pfc-exp-cli` is packaged alongside the GUI.");
        return;
    }
#ifdef Q_OS_WIN
    {
        const QDir exeDir(QFileInfo(program).absolutePath());
        const QStringList fftwDlls = exeDir.entryList(QStringList() << "*fftw*.dll", QDir::Files);
        if (fftwDlls.isEmpty()) {
            appendLog("=== Warning: no *fftw*.dll found next to pfc-exp-cli; Windows may fail to start (0xC0000135) ===");
        }
    }
#endif
#ifdef Q_OS_MACOS
    {
        const QDir exeDir(QCoreApplication::applicationDirPath());
        const QDir frameworksDir(exeDir.filePath("../Frameworks"));
        if (frameworksDir.exists()) {
            const QStringList fftwDylibs = frameworksDir.entryList(QStringList() << "libfftw3*.dylib", QDir::Files);
            if (fftwDylibs.isEmpty()) {
                appendLog("=== Warning: no libfftw3*.dylib found in Contents/Frameworks; macOS may crash at start if FFTW is missing ===");
            }
        }
    }
#endif

    if (process_) {
        process_->deleteLater();
        process_ = nullptr;
    }

    process_ = new QProcess(this);
    process_->setProgram(program);
    process_->setArguments(buildArgs(job.params, job.outDir));
    process_->setWorkingDirectory(job.outDir);
    process_->setProcessChannelMode(QProcess::MergedChannels);

    connect(process_, &QProcess::readyRead, this, &MainWindow::onProcessReadyRead);
    connect(process_, &QProcess::finished, this, &MainWindow::onProcessFinished);
    connect(process_, &QProcess::errorOccurred, this, &MainWindow::onProcessError);

    currentOutputDir_ = job.outDir;
    appendLog(QString("=== %1: %2 ===").arg(formatRunLabel(job.index, job.total), job.outDir));
    appendLog(program + " " + process_->arguments().join(' '));

    process_->start();
}

void MainWindow::startSingleRun() {
    if (process_) {
        QMessageBox::warning(this, "Busy", "A run is already in progress.");
        return;
    }

    QString err;
    if (!ensureDir(singleBaseOutDir_->text(), &err)) {
        QMessageBox::critical(this, "Output error", err);
        return;
    }

    const QString runName = singleRunName_->text().trimmed().isEmpty() ? makeTimestampedDirName("single") : singleRunName_->text().trimmed();
    const QString outDir = QDir(singleBaseOutDir_->text()).filePath(runName);

    RunJob job;
    job.mode = "single";
    job.index = 1;
    job.total = 1;
    job.params = singleParams();
    job.outDir = outDir;
    launchJob(job);
}

void MainWindow::startBatchRun() {
    if (process_) {
        QMessageBox::warning(this, "Busy", "A run is already in progress.");
        return;
    }
    const QString baseDir = batchBaseOutDir_->text().trimmed();
    const QString batchName = batchName_->text().trimmed().isEmpty() ? makeTimestampedDirName("batch") : batchName_->text().trimmed();

    QString err;
    jobQueue_ = buildBatchJobs(baseDir, batchName, &err);
    if (!err.isEmpty()) {
        QMessageBox::critical(this, "Batch config error", err);
        return;
    }
    if (jobQueue_.isEmpty()) {
        QMessageBox::warning(this, "Batch", "No runs to execute.");
        return;
    }
    currentJobIndex_ = 0;
    batchProgress_->setValue(0);
    launchJob(jobQueue_.at(currentJobIndex_));
}

void MainWindow::stopRun() {
    jobQueue_.clear();
    currentJobIndex_ = -1;
    if (!process_) return;
    appendLog("=== Stop requested ===");
    process_->kill();
}

void MainWindow::browseSingleOutDir() {
    const QString dir = QFileDialog::getExistingDirectory(this, "Select output directory", singleBaseOutDir_->text());
    if (!dir.isEmpty()) singleBaseOutDir_->setText(dir);
}

void MainWindow::browseBatchOutDir() {
    const QString dir = QFileDialog::getExistingDirectory(this, "Select output directory", batchBaseOutDir_->text());
    if (!dir.isEmpty()) batchBaseOutDir_->setText(dir);
}

void MainWindow::openCurrentOutputDir() {
    const QString dir = currentOutputDir_.isEmpty() ? singleBaseOutDir_->text() : currentOutputDir_;
    if (dir.trimmed().isEmpty()) return;
    QDesktopServices::openUrl(QUrl::fromLocalFile(dir));
}

void MainWindow::updateBatchPreview() {
    QString err;
    const QVector<double> u0Values = sweepU0Enable_->isChecked()
        ? makeDoubleRange(sweepU0Start_->value(), sweepU0End_->value(), sweepU0Step_->value(), &err)
        : QVector<double>{sweepU0Single_->value()};
    if (!err.isEmpty()) {
        batchPreview_->setText("Batch preview: u0: " + err);
        return;
    }
    const QVector<double> con0Values = sweepCon0Enable_->isChecked()
        ? makeDoubleRange(sweepCon0Start_->value(), sweepCon0End_->value(), sweepCon0Step_->value(), &err)
        : QVector<double>{sweepCon0Single_->value()};
    if (!err.isEmpty()) {
        batchPreview_->setText("Batch preview: con0: " + err);
        return;
    }
    const QVector<int> stepsValues = sweepStepsEnable_->isChecked()
        ? makeIntRange(sweepStepsStart_->value(), sweepStepsEnd_->value(), sweepStepsStep_->value(), &err)
        : QVector<int>{sweepStepsSingle_->value()};
    if (!err.isEmpty()) {
        batchPreview_->setText("Batch preview: steps: " + err);
        return;
    }

    const long long total =
        static_cast<long long>(u0Values.size()) * static_cast<long long>(con0Values.size()) * static_cast<long long>(stepsValues.size());
    batchPreview_->setText(QString("Batch preview: %1 run(s) will be executed.").arg(total));
}

void MainWindow::onProcessReadyRead() {
    if (!process_) return;
    const QByteArray data = process_->readAll();
    if (!data.isEmpty()) appendLog(QString::fromLocal8Bit(data));
}

void MainWindow::onProcessFinished(int exitCode, QProcess::ExitStatus exitStatus) {
    const bool ok = (exitStatus == QProcess::NormalExit && exitCode == 0);
    const QString statusStr = (exitStatus == QProcess::NormalExit) ? "NormalExit" : "CrashExit";
    const quint32 code = static_cast<quint32>(exitCode);
    const QString hexCode = QString("0x%1").arg(code, 8, 16, QChar('0')).toUpper();
    appendLog(QString("=== Finished: exitStatus=%1 exitCode=%2 (%3, %4) ===").arg(statusStr).arg(exitCode).arg(hexCode).arg(ok ? "OK" : "FAILED"));
#ifdef Q_OS_WIN
    if (!ok && exitStatus == QProcess::CrashExit && code == 0xC0000135u) {
        appendLog("=== Hint: 0xC0000135 usually means a required DLL was not found. Ensure FFTW DLLs are bundled next to pfc-exp-cli.exe ===");
    }
#endif
#ifdef Q_OS_MACOS
    if (!ok && exitStatus == QProcess::CrashExit && code == 9u) {
        appendLog("=== Hint: SIGKILL(9) on macOS can be Gatekeeper/quarantine or a missing dylib. If this is a downloaded app, try: xattr -dr com.apple.quarantine /Applications/PFC-Exp-GUI.app ===");
    }
#endif

    process_->deleteLater();
    process_ = nullptr;

    if (jobQueue_.isEmpty() || currentJobIndex_ < 0) {
        batchProgress_->setValue(0);
        currentJobIndex_ = -1;
        return;
    }

    const int total = jobQueue_.size();
    const int done = currentJobIndex_ + 1;
    batchProgress_->setValue(static_cast<int>(100.0 * done / total));

    if (!ok) {
        appendLog("=== Batch aborted due to failure ===");
        jobQueue_.clear();
        currentJobIndex_ = -1;
        return;
    }

    ++currentJobIndex_;
    if (currentJobIndex_ >= total) {
        appendLog("=== Batch complete ===");
        jobQueue_.clear();
        currentJobIndex_ = -1;
        batchProgress_->setValue(100);
        return;
    }

    launchJob(jobQueue_.at(currentJobIndex_));
}

void MainWindow::onProcessError(QProcess::ProcessError error) {
    if (!process_) {
        appendLog(QString("=== Process error: %1 (process=null) ===").arg(static_cast<int>(error)));
        return;
    }

    QString errorName;
    switch (error) {
        case QProcess::FailedToStart: errorName = "FailedToStart"; break;
        case QProcess::Crashed: errorName = "Crashed"; break;
        case QProcess::Timedout: errorName = "Timedout"; break;
        case QProcess::WriteError: errorName = "WriteError"; break;
        case QProcess::ReadError: errorName = "ReadError"; break;
        case QProcess::UnknownError: errorName = "UnknownError"; break;
    }

    appendLog(QString("=== Process error: %1 (%2) ===").arg(static_cast<int>(error)).arg(errorName));
    appendLog("=== errorString: " + process_->errorString() + " ===");
    appendLog("=== program: " + process_->program() + " ===");
    appendLog("=== workingDir: " + process_->workingDirectory() + " ===");
}

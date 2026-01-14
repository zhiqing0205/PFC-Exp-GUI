#pragma once

#include <QObject>
#include <QFile>
#include <QProcess>
#include <QQueue>
#include <QString>
#include <QVariantMap>

#include "RunQueueModel.h"

class ExperimentController final : public QObject {
  Q_OBJECT

  Q_PROPERTY(RunQueueModel* runModel READ runModel CONSTANT)
  Q_PROPERTY(bool running READ running NOTIFY runningChanged)
  Q_PROPERTY(int totalRuns READ totalRuns NOTIFY totalsChanged)
  Q_PROPERTY(int finishedRuns READ finishedRuns NOTIFY totalsChanged)
  Q_PROPERTY(int selectedRun READ selectedRun WRITE setSelectedRun NOTIFY selectedRunChanged)
  Q_PROPERTY(QString selectedLog READ selectedLog NOTIFY selectedLogChanged)
  Q_PROPERTY(QString lastError READ lastError NOTIFY lastErrorChanged)

public:
  explicit ExperimentController(QObject* parent = nullptr);

  RunQueueModel* runModel();

  bool running() const;
  int totalRuns() const;
  int finishedRuns() const;

  int selectedRun() const;
  void setSelectedRun(int row);

  QString selectedLog() const;
  QString lastError() const;

  Q_INVOKABLE void startRuns(const QVariantMap& config);
  Q_INVOKABLE void stop();
  Q_INVOKABLE void clear();

signals:
  void runningChanged();
  void totalsChanged();
  void selectedRunChanged();
  void selectedLogChanged();
  void lastErrorChanged();

private:
  struct RunSpec {
    int row = -1;
    QString name;
    QString outDir;
    QString program;
    QStringList args;
    QString commandString;
    QString guiLogPath;
  };

  void setRunning(bool v);
  void setLastError(const QString& err);
  void updateSelectedLog();

  void startNext();
  void handleProcessOutput();
  void handleProcessFinished(int exitCode, QProcess::ExitStatus exitStatus);
  void handleProcessError(QProcess::ProcessError error);

  static QString defaultSimPath();
  static QString normalizedPath(QString p);
  static QString safeToken(QString s);
  static QVector<double> buildDoubleValues(bool enabled, double single, double start, double end, double step);
  static QVector<int> buildIntValues(bool enabled, int single, int start, int end, int step);

  RunQueueModel model_;
  QQueue<RunSpec> queue_;
  QProcess process_;
  QFile currentGuiLog_;

  bool running_ = false;
  int totalRuns_ = 0;
  int finishedRuns_ = 0;
  int selectedRun_ = -1;
  QString selectedLog_;
  QString lastError_;
  int currentRow_ = -1;
};

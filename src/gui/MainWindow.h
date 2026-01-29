#pragma once

#include <QMainWindow>
#include <QProcess>
#include <QVector>
#include <QSize>

class QCheckBox;
class QDoubleSpinBox;
class QComboBox;
class QGraphicsPixmapItem;
class QGraphicsScene;
class QGraphicsView;
class QLabel;
class QLineEdit;
class QListWidget;
class QPlainTextEdit;
class QProgressBar;
class QSpinBox;
class QStackedWidget;
class QTabWidget;
class QTableWidget;

class PlotWidget;

struct RunParams {
    double u0 = 0.05;
    double con0 = 0.2;
    double sig = 0.05;
    double dt = 0.05;
    double dx = 0.125;
    int steps = 200;
    int mod = 25;
    int seed = 20200604;
};

struct RunJob {
    RunParams params;
    QString outDir;
    int index = 0;
    int total = 1;
    QString mode;
};

class MainWindow final : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

private slots:
    void startSingleRun();
    void startBatchRun();
    void stopRun();
    void browseSingleOutDir();
    void browseBatchOutDir();
    void openCurrentOutputDir();
    void updateBatchPreview();

    void onProcessReadyRead();
    void onProcessFinished(int exitCode, QProcess::ExitStatus exitStatus);
    void onProcessError(QProcess::ProcessError error);

private:
    QString cliPath() const;
    QString defaultBaseOutputDir() const;
    bool ensureDir(const QString& dirPath, QString* errorOut) const;
    QString makeTimestampedDirName(const QString& prefix) const;
    QString formatRunLabel(int index, int total) const;

    RunParams singleParams() const;
    QVector<RunJob> buildBatchJobs(const QString& baseDir, const QString& batchName, QString* errorOut) const;
    QStringList buildArgs(const RunParams& params, const QString& outDir) const;
    bool writeParamsJson(const RunJob& job, QString* errorOut) const;

    void launchJob(const RunJob& job);
    void appendLog(const QString& text);
    void setRunningUi(bool running);

    void setResultsDir(const QString& dirPath);
    void refreshResultsFileList();
    void loadResultsFile(const QString& filePath);
    void applyResultsYColumn(int columnIndex);
    void renderResultsImage();
    void fitResultsImage();
    void exportResultsImage();

    QProcess* process_ = nullptr;
    QVector<RunJob> jobQueue_;
    int currentJobIndex_ = -1;
    QString currentOutputDir_;
    QString processOutputBuffer_;
    QString currentJobMode_;
    int currentJobTotalSteps_ = 0;

    QPlainTextEdit* log_ = nullptr;
    QProgressBar* singleStepProgress_ = nullptr;
    QProgressBar* batchStepProgress_ = nullptr;

    // Single-run widgets
    QDoubleSpinBox* singleU0_ = nullptr;
    QDoubleSpinBox* singleCon0_ = nullptr;
    QDoubleSpinBox* singleSig_ = nullptr;
    QDoubleSpinBox* singleDt_ = nullptr;
    QDoubleSpinBox* singleDx_ = nullptr;
    QSpinBox* singleSteps_ = nullptr;
    QSpinBox* singleMod_ = nullptr;
    QSpinBox* singleSeed_ = nullptr;
    QLineEdit* singleBaseOutDir_ = nullptr;
    QLineEdit* singleRunName_ = nullptr;

    // Batch-run widgets
    struct BatchParamWidgets {
        QString key;
        QString label;
        bool isInt = false;
        int decimals = 6;

        QComboBox* mode = nullptr;         // Fixed / Range / List
        QStackedWidget* stack = nullptr;   // page per mode

        // Fixed
        QSpinBox* fixedInt = nullptr;
        QDoubleSpinBox* fixedDouble = nullptr;

        // Range
        QSpinBox* rangeStartInt = nullptr;
        QSpinBox* rangeEndInt = nullptr;
        QSpinBox* rangeStepInt = nullptr;
        QDoubleSpinBox* rangeStartDouble = nullptr;
        QDoubleSpinBox* rangeEndDouble = nullptr;
        QDoubleSpinBox* rangeStepDouble = nullptr;
        QLineEdit* rangePreview = nullptr;  // read-only, comma-separated list

        // List
        QLineEdit* listEdit = nullptr;      // comma-separated values
    };

    QVector<BatchParamWidgets> batchParams_;

    QLineEdit* batchBaseOutDir_ = nullptr;
    QLineEdit* batchName_ = nullptr;
    QLabel* batchPreview_ = nullptr;
    QProgressBar* batchProgress_ = nullptr;

    // Results/visualization widgets
    QLineEdit* resultsDir_ = nullptr;
    QListWidget* resultsFiles_ = nullptr;
    QTabWidget* resultsViewTabs_ = nullptr;
    QWidget* resultsImagePage_ = nullptr;
    QComboBox* resultsYColumn_ = nullptr;
    QLabel* resultsInfo_ = nullptr;
    PlotWidget* resultsPlot_ = nullptr;
    QTableWidget* resultsTable_ = nullptr;
    QComboBox* resultsImageValue_ = nullptr;
    QSpinBox* resultsImagePointSize_ = nullptr;
    QLabel* resultsImageInfo_ = nullptr;
    QGraphicsView* resultsImageView_ = nullptr;
    QGraphicsScene* resultsImageScene_ = nullptr;
    QGraphicsPixmapItem* resultsImageItem_ = nullptr;
    QString resultsCurrentFile_;
    bool resultsImageEnabled_ = false;
    QVector<QVector<double>> resultsColumns_;

    // Heatmap rendering buffers (reusable to avoid frequent allocations)
    QSize resultsImageBufferSize_;
    QVector<float> resultsImageSum_;
    QVector<float> resultsImageWsum_;
    QVector<float> resultsImageWmax_;
};

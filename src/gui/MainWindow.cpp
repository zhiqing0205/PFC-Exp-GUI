#include "MainWindow.h"

#include "PlotWidget.h"

#include <QCheckBox>
#include <QComboBox>
#include <QCoreApplication>
#include <QDateTime>
#include <QDesktopServices>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QFrame>
#include <QGridLayout>
#include <QGroupBox>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QHBoxLayout>
#include <QHash>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPainter>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QRegularExpression>
#include <QTextStream>
#include <QScrollArea>
#include <QSplitter>
#include <QSpinBox>
#include <QStandardPaths>
#include <QStyle>
#include <QTabBar>
#include <QTabWidget>
#include <QToolButton>
#include <QUrl>
#include <QVBoxLayout>
#include <QFontDatabase>
#include <QTableWidget>
#include <QHeaderView>
#include <QWheelEvent>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>

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

class ZoomableGraphicsView final : public QGraphicsView {
public:
    using QGraphicsView::QGraphicsView;

protected:
    void wheelEvent(QWheelEvent* event) override {
        if (event->modifiers().testFlag(Qt::ControlModifier)) {
            constexpr double base = 1.0015;
            const double factor = std::pow(base, static_cast<double>(event->angleDelta().y()));
            scale(factor, factor);
            event->accept();
            return;
        }
        QGraphicsView::wheelEvent(event);
    }
};

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("PFC-Exp-GUI");
    resize(1100, 780);

    auto* tabs = new QTabWidget;
    tabs->setDocumentMode(true);
    tabs->setUsesScrollButtons(true);
    tabs->tabBar()->setDrawBase(false);
    tabs->tabBar()->setExpanding(false);

    auto* singleTab = new QWidget;
    auto* batchTab = new QWidget;
    auto* resultsTab = new QWidget;
    auto* manufacturingTab = new QWidget;
    auto* transformationTab = new QWidget;
    auto* mechanicsTab = new QWidget;

    auto* singleScroll = new QScrollArea;
    singleScroll->setFrameShape(QFrame::NoFrame);
    singleScroll->setWidgetResizable(true);
    singleScroll->setWidget(singleTab);

    auto* batchScroll = new QScrollArea;
    batchScroll->setFrameShape(QFrame::NoFrame);
    batchScroll->setWidgetResizable(true);
    batchScroll->setWidget(batchTab);

    auto* resultsScroll = new QScrollArea;
    resultsScroll->setFrameShape(QFrame::NoFrame);
    resultsScroll->setWidgetResizable(true);
    resultsScroll->setWidget(resultsTab);

    auto* manufacturingScroll = new QScrollArea;
    manufacturingScroll->setFrameShape(QFrame::NoFrame);
    manufacturingScroll->setWidgetResizable(true);
    manufacturingScroll->setWidget(manufacturingTab);

    auto* transformationScroll = new QScrollArea;
    transformationScroll->setFrameShape(QFrame::NoFrame);
    transformationScroll->setWidgetResizable(true);
    transformationScroll->setWidget(transformationTab);

    auto* mechanicsScroll = new QScrollArea;
    mechanicsScroll->setFrameShape(QFrame::NoFrame);
    mechanicsScroll->setWidgetResizable(true);
    mechanicsScroll->setWidget(mechanicsTab);

    tabs->addTab(singleScroll, "Single Run");
    tabs->addTab(batchScroll, "Batch Sweep");
    tabs->addTab(resultsScroll, "TXT Visualizer");
    tabs->addTab(manufacturingScroll, "Manufacturing");
    tabs->addTab(transformationScroll, "Transformation");
    tabs->addTab(mechanicsScroll, "Stress–Strain");

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
    auto* splitter = new QSplitter(Qt::Vertical);
    splitter->setHandleWidth(10);
    splitter->setChildrenCollapsible(false);
    splitter->addWidget(tabs);
    splitter->addWidget(logBox);
    splitter->setStretchFactor(0, 3);
    splitter->setStretchFactor(1, 1);
    splitter->setSizes(QList<int>() << 560 << 220);

    centralLayout->addWidget(splitter, 1);
    central->setLayout(centralLayout);
    setCentralWidget(central);

    // ---------- Single tab ----------
    auto* singleLayout = new QVBoxLayout;
    singleLayout->setSpacing(12);
    singleTab->setLayout(singleLayout);

    auto* singleParamsBox = new QGroupBox("Parameters");
    auto* singleParamsGrid = new QGridLayout;
    singleParamsGrid->setVerticalSpacing(8);
    singleParamsGrid->setHorizontalSpacing(12);
    singleParamsGrid->setColumnStretch(1, 1);
    singleParamsGrid->setColumnStretch(3, 1);

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

    auto addSingleParam = [&](int row, int colPair, const QString& labelText, QWidget* editor) {
        auto* label = new QLabel(labelText);
        label->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        singleParamsGrid->addWidget(label, row, colPair * 2);
        singleParamsGrid->addWidget(editor, row, colPair * 2 + 1);
    };

    addSingleParam(0, 0, "u0 (density)", singleU0_);
    addSingleParam(0, 1, "con0 (concentration)", singleCon0_);
    addSingleParam(1, 0, "sig", singleSig_);
    addSingleParam(1, 1, "seed", singleSeed_);
    addSingleParam(2, 0, "dt", singleDt_);
    addSingleParam(2, 1, "dx", singleDx_);
    addSingleParam(3, 0, "steps", singleSteps_);
    addSingleParam(3, 1, "mod", singleMod_);
    addSingleParam(4, 0, "grainx", singleGrainX_);
    addSingleParam(4, 1, "grainy", singleGrainY_);
    addSingleParam(5, 0, "grainz", singleGrainZ_);
    addSingleParam(5, 1, "axx (noise)", singleAxx_);

    singleParamsBox->setLayout(singleParamsGrid);

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

    auto* batchHint = new QLabel(
        "Batch Sweep: enable Sweep for u0 / con0 / steps, then fill Start/End/Step.\n"
        "Disabled fields are ignored.");
    batchHint->setProperty("hint", true);
    batchHint->setWordWrap(true);

    auto* batchSweepBox = new QGroupBox("Sweep (Cartesian product)");
    auto* sweepGrid = new QGridLayout;
    sweepGrid->setVerticalSpacing(8);
    sweepGrid->setHorizontalSpacing(10);
    sweepGrid->setColumnStretch(2, 1);
    sweepGrid->setColumnStretch(3, 1);
    sweepGrid->setColumnStretch(4, 1);
    sweepGrid->setColumnStretch(5, 1);

    auto header = [](const QString& text) {
        auto* l = new QLabel(text);
        QFont f = l->font();
        f.setBold(true);
        l->setFont(f);
        l->setAlignment(Qt::AlignCenter);
        return l;
    };

    sweepGrid->addWidget(header("Param"), 0, 0);
    sweepGrid->addWidget(header("Sweep"), 0, 1);
    sweepGrid->addWidget(header("Single"), 0, 2);
    sweepGrid->addWidget(header("Start"), 0, 3);
    sweepGrid->addWidget(header("End"), 0, 4);
    sweepGrid->addWidget(header("Step"), 0, 5);

    sweepU0Enable_ = new QCheckBox;
    sweepU0Enable_->setToolTip("Enable sweep for u0");
    sweepU0Single_ = makeDoubleSpin(-2.0, 2.0, 0.05, 6, 0.01);
    sweepU0Start_ = makeDoubleSpin(-2.0, 2.0, 0.0, 6, 0.01);
    sweepU0End_ = makeDoubleSpin(-2.0, 2.0, 0.2, 6, 0.01);
    sweepU0Step_ = makeDoubleSpin(1e-9, 10.0, 0.05, 6, 0.01);

    sweepCon0Enable_ = new QCheckBox;
    sweepCon0Enable_->setToolTip("Enable sweep for con0");
    sweepCon0Single_ = makeDoubleSpin(0.0, 1.0, 0.2, 6, 0.01);
    sweepCon0Start_ = makeDoubleSpin(0.0, 1.0, 0.1, 6, 0.01);
    sweepCon0End_ = makeDoubleSpin(0.0, 1.0, 0.3, 6, 0.01);
    sweepCon0Step_ = makeDoubleSpin(1e-9, 1.0, 0.05, 6, 0.01);

    sweepStepsEnable_ = new QCheckBox;
    sweepStepsEnable_->setToolTip("Enable sweep for steps");
    sweepStepsSingle_ = makeIntSpin(1, 100000000, 200, 10);
    sweepStepsStart_ = makeIntSpin(1, 100000000, 100, 10);
    sweepStepsEnd_ = makeIntSpin(1, 100000000, 300, 10);
    sweepStepsStep_ = makeIntSpin(1, 100000000, 100, 10);

    auto addSweepRow = [&](int row, const QString& name, QCheckBox* enable, QWidget* single, QWidget* start, QWidget* end, QWidget* step) {
        auto* label = new QLabel(name);
        label->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        sweepGrid->addWidget(label, row, 0);
        sweepGrid->addWidget(enable, row, 1, Qt::AlignCenter);
        sweepGrid->addWidget(single, row, 2);
        sweepGrid->addWidget(start, row, 3);
        sweepGrid->addWidget(end, row, 4);
        sweepGrid->addWidget(step, row, 5);
    };

    addSweepRow(1, "u0", sweepU0Enable_, sweepU0Single_, sweepU0Start_, sweepU0End_, sweepU0Step_);
    addSweepRow(2, "con0", sweepCon0Enable_, sweepCon0Single_, sweepCon0Start_, sweepCon0End_, sweepCon0Step_);
    addSweepRow(3, "steps", sweepStepsEnable_, sweepStepsSingle_, sweepStepsStart_, sweepStepsEnd_, sweepStepsStep_);

    batchSweepBox->setLayout(sweepGrid);

    auto* batchFixedBox = new QGroupBox("Fixed Parameters");
    auto* batchFixedGrid = new QGridLayout;
    batchFixedGrid->setVerticalSpacing(8);
    batchFixedGrid->setHorizontalSpacing(12);
    batchFixedGrid->setColumnStretch(1, 1);
    batchFixedGrid->setColumnStretch(3, 1);
    batchSig_ = makeDoubleSpin(0.0, 2.0, 0.05, 6, 0.01);
    batchDt_ = makeDoubleSpin(1e-6, 1.0, 0.05, 6, 0.01);
    batchDx_ = makeDoubleSpin(1e-6, 10.0, 0.125, 6, 0.01);
    batchMod_ = makeIntSpin(1, 100000000, 25, 1);
    batchSeed_ = makeIntSpin(0, 2147483647, 20200604, 1);
    batchGrainX_ = makeIntSpin(1, 4096, 16, 1);
    batchGrainY_ = makeIntSpin(1, 4096, 16, 1);
    batchGrainZ_ = makeIntSpin(1, 4096, 1, 1);
    batchAxx_ = makeDoubleSpin(0.0, 1000.0, 0.0, 6, 0.1);

    auto addBatchFixed = [&](int row, int colPair, const QString& labelText, QWidget* editor) {
        auto* label = new QLabel(labelText);
        label->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        batchFixedGrid->addWidget(label, row, colPair * 2);
        batchFixedGrid->addWidget(editor, row, colPair * 2 + 1);
    };

    addBatchFixed(0, 0, "sig", batchSig_);
    addBatchFixed(0, 1, "seed", batchSeed_);
    addBatchFixed(1, 0, "dt", batchDt_);
    addBatchFixed(1, 1, "dx", batchDx_);
    addBatchFixed(2, 0, "mod", batchMod_);
    addBatchFixed(2, 1, "axx (noise)", batchAxx_);
    addBatchFixed(3, 0, "grainx", batchGrainX_);
    addBatchFixed(3, 1, "grainy", batchGrainY_);
    addBatchFixed(4, 0, "grainz", batchGrainZ_);

    batchFixedBox->setLayout(batchFixedGrid);

    auto updateSweepEnabled = [this]() {
        const bool u0On = sweepU0Enable_->isChecked();
        sweepU0Single_->setEnabled(!u0On);
        sweepU0Start_->setEnabled(u0On);
        sweepU0End_->setEnabled(u0On);
        sweepU0Step_->setEnabled(u0On);

        const bool con0On = sweepCon0Enable_->isChecked();
        sweepCon0Single_->setEnabled(!con0On);
        sweepCon0Start_->setEnabled(con0On);
        sweepCon0End_->setEnabled(con0On);
        sweepCon0Step_->setEnabled(con0On);

        const bool stepsOn = sweepStepsEnable_->isChecked();
        sweepStepsSingle_->setEnabled(!stepsOn);
        sweepStepsStart_->setEnabled(stepsOn);
        sweepStepsEnd_->setEnabled(stepsOn);
        sweepStepsStep_->setEnabled(stepsOn);
    };

    connect(sweepU0Enable_, &QCheckBox::toggled, this, [=](bool) {
        updateSweepEnabled();
        updateBatchPreview();
    });
    connect(sweepCon0Enable_, &QCheckBox::toggled, this, [=](bool) {
        updateSweepEnabled();
        updateBatchPreview();
    });
    connect(sweepStepsEnable_, &QCheckBox::toggled, this, [=](bool) {
        updateSweepEnabled();
        updateBatchPreview();
    });

    updateSweepEnabled();

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

    batchLayout->addWidget(batchHint);
    batchLayout->addWidget(batchSweepBox);
    batchLayout->addWidget(batchPreview_);
    batchLayout->addWidget(batchFixedBox);
    batchLayout->addWidget(batchOutBox);
    batchLayout->addWidget(batchProgress_);
    batchLayout->addWidget(batchButtons);
    batchLayout->addStretch(1);

    // ---------- Results tab ----------
    auto* resultsLayout = new QVBoxLayout;
    resultsLayout->setSpacing(12);
    resultsTab->setLayout(resultsLayout);

    auto* resultsHint = new QLabel(
        "Select an output run directory to visualize numeric data in *.txt files.\n"
        "Tips:\n"
        "• energy.txt / phimax.txt: use Plot (x=step)\n"
        "• Phimax_*.txt: use Image to render a 2D snapshot (Ctrl+wheel to zoom)");
    resultsHint->setProperty("hint", true);
    resultsHint->setWordWrap(true);
    resultsLayout->addWidget(resultsHint);

    auto* resultsDirBox = new QGroupBox("Run directory");
    auto* resultsDirRow = new QHBoxLayout;
    resultsDirRow->setContentsMargins(0, 0, 0, 0);
    resultsDir_ = new QLineEdit(defaultBaseOutputDir());
    resultsDir_->setPlaceholderText("Select a run directory containing *.txt outputs");
    auto* resultsBrowse = new QToolButton;
    resultsBrowse->setText("Browse…");
    resultsBrowse->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    resultsBrowse->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    auto* resultsRefresh = new QToolButton;
    resultsRefresh->setToolButtonStyle(Qt::ToolButtonIconOnly);
    resultsRefresh->setIcon(style()->standardIcon(QStyle::SP_BrowserReload));
    resultsRefresh->setToolTip("Refresh file list");
    auto* resultsUseCurrent = new QPushButton("Use Current Output");
    resultsUseCurrent->setIcon(style()->standardIcon(QStyle::SP_DialogApplyButton));

    resultsDirRow->addWidget(resultsDir_, 1);
    resultsDirRow->addWidget(resultsBrowse);
    resultsDirRow->addWidget(resultsRefresh);
    resultsDirRow->addWidget(resultsUseCurrent);
    resultsDirBox->setLayout(resultsDirRow);
    resultsLayout->addWidget(resultsDirBox);

    auto* resultsSplit = new QSplitter(Qt::Horizontal);
    resultsSplit->setHandleWidth(10);
    resultsSplit->setChildrenCollapsible(false);

    auto* filesBox = new QGroupBox("TXT files");
    auto* filesLayout = new QVBoxLayout;
    filesLayout->setContentsMargins(0, 0, 0, 0);
    resultsFiles_ = new QListWidget;
    resultsFiles_->setMinimumWidth(220);
    filesLayout->addWidget(resultsFiles_);
    filesBox->setLayout(filesLayout);

    resultsViewTabs_ = new QTabWidget;
    resultsViewTabs_->setDocumentMode(true);

    auto* plotPage = new QWidget;
    auto* plotLayout = new QVBoxLayout;
    plotLayout->setContentsMargins(0, 0, 0, 0);
    plotLayout->setSpacing(10);

    auto* plotTopRow = new QWidget;
    auto* plotTopLayout = new QHBoxLayout;
    plotTopLayout->setContentsMargins(0, 0, 0, 0);
    plotTopLayout->addWidget(new QLabel("Y column"));
    resultsYColumn_ = new QComboBox;
    resultsYColumn_->setMinimumWidth(180);
    plotTopLayout->addWidget(resultsYColumn_);
    plotTopLayout->addStretch(1);
    plotTopRow->setLayout(plotTopLayout);

    resultsPlot_ = new PlotWidget;
    resultsInfo_ = new QLabel("Pick a file to visualize.");
    resultsInfo_->setProperty("hint", true);
    resultsInfo_->setWordWrap(true);

    plotLayout->addWidget(plotTopRow);
    plotLayout->addWidget(resultsPlot_, 1);
    plotLayout->addWidget(resultsInfo_);
    plotPage->setLayout(plotLayout);

    auto* tablePage = new QWidget;
    auto* tableLayout = new QVBoxLayout;
    tableLayout->setContentsMargins(0, 0, 0, 0);
    resultsTable_ = new QTableWidget;
    resultsTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    resultsTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
    resultsTable_->setAlternatingRowColors(true);
    resultsTable_->horizontalHeader()->setStretchLastSection(true);
    tableLayout->addWidget(resultsTable_);
    tablePage->setLayout(tableLayout);

    auto* imagePage = new QWidget;
    auto* imageLayout = new QVBoxLayout;
    imageLayout->setContentsMargins(0, 0, 0, 0);
    imageLayout->setSpacing(10);

    auto* imageTopRow = new QWidget;
    auto* imageTopLayout = new QHBoxLayout;
    imageTopLayout->setContentsMargins(0, 0, 0, 0);
    imageTopLayout->addWidget(new QLabel("Value"));
    resultsImageValue_ = new QComboBox;
    resultsImageValue_->setMinimumWidth(180);
    resultsImageValue_->setEnabled(false);
    imageTopLayout->addWidget(resultsImageValue_);
    imageTopLayout->addWidget(new QLabel("Smooth σ (px)"));
    resultsImagePointSize_ = makeIntSpin(1, 64, 16, 1);
    resultsImagePointSize_->setToolTip("Gaussian smoothing radius (σ). Higher values create smoother heatmaps.");
    resultsImagePointSize_->setEnabled(false);
    imageTopLayout->addWidget(resultsImagePointSize_);
    imageTopLayout->addStretch(1);

    auto* imageFit = new QToolButton;
    imageFit->setText("Fit");
    imageFit->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    imageFit->setIcon(style()->standardIcon(QStyle::SP_DesktopIcon));
    imageTopLayout->addWidget(imageFit);

    auto* imageExport = new QToolButton;
    imageExport->setText("Export PNG…");
    imageExport->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    imageExport->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
    imageTopLayout->addWidget(imageExport);
    imageTopRow->setLayout(imageTopLayout);

    resultsImageInfo_ = new QLabel("Select a Phimax_*.txt snapshot file to render.");
    resultsImageInfo_->setProperty("hint", true);
    resultsImageInfo_->setWordWrap(true);

    resultsImageScene_ = new QGraphicsScene(this);
    resultsImageView_ = new ZoomableGraphicsView;
    resultsImageView_->setScene(resultsImageScene_);
    resultsImageView_->setFrameShape(QFrame::NoFrame);
    resultsImageView_->setBackgroundBrush(QColor(248, 248, 250));
    resultsImageView_->setRenderHint(QPainter::Antialiasing, true);
    resultsImageView_->setDragMode(QGraphicsView::ScrollHandDrag);
    resultsImageView_->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    resultsImageView_->setResizeAnchor(QGraphicsView::AnchorUnderMouse);
    resultsImageItem_ = resultsImageScene_->addPixmap(QPixmap());

    imageLayout->addWidget(imageTopRow);
    imageLayout->addWidget(resultsImageView_, 1);
    imageLayout->addWidget(resultsImageInfo_);
    imagePage->setLayout(imageLayout);
    resultsImagePage_ = imagePage;

    resultsViewTabs_->addTab(plotPage, "Plot");
    resultsViewTabs_->addTab(tablePage, "Table");
    resultsViewTabs_->addTab(imagePage, "Image");

    resultsSplit->addWidget(filesBox);
    resultsSplit->addWidget(resultsViewTabs_);
    resultsSplit->setStretchFactor(0, 1);
    resultsSplit->setStretchFactor(1, 3);
    resultsSplit->setSizes(QList<int>() << 260 << 780);

    resultsLayout->addWidget(resultsSplit, 1);

    connect(resultsBrowse, &QToolButton::clicked, this, [this]() {
        const QString base = resultsDir_ ? resultsDir_->text() : defaultBaseOutputDir();
        const QString dir = QFileDialog::getExistingDirectory(this, "Select run directory", base);
        if (!dir.isEmpty()) setResultsDir(dir);
    });
    connect(resultsRefresh, &QToolButton::clicked, this, [this]() { refreshResultsFileList(); });
    connect(resultsUseCurrent, &QPushButton::clicked, this, [this]() {
        if (!currentOutputDir_.isEmpty()) setResultsDir(currentOutputDir_);
    });
    connect(resultsDir_, &QLineEdit::editingFinished, this, [this]() { refreshResultsFileList(); });
    connect(resultsFiles_, &QListWidget::currentTextChanged, this, [this](const QString& fileName) {
        if (fileName.trimmed().isEmpty()) return;
        const QString dir = resultsDir_ ? resultsDir_->text().trimmed() : QString();
        if (dir.isEmpty()) return;
        loadResultsFile(QDir(dir).filePath(fileName));
    });
    connect(resultsYColumn_, &QComboBox::currentIndexChanged, this, [this](int) {
        if (!resultsYColumn_) return;
        const int col = resultsYColumn_->currentData().toInt();
        applyResultsYColumn(col);
    });

    connect(resultsImageValue_, &QComboBox::currentIndexChanged, this, [this](int) { renderResultsImage(); });
    connect(resultsImagePointSize_, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int) { renderResultsImage(); });
    connect(imageFit, &QToolButton::clicked, this, [this]() { fitResultsImage(); });
    connect(imageExport, &QToolButton::clicked, this, [this]() { exportResultsImage(); });

    // ---------- Roadmap tabs (placeholders) ----------
    auto setupRoadmapTab = [](QWidget* tab, const QString& title, const QString& desc, const QStringList& bullets) {
        auto* layout = new QVBoxLayout;
        layout->setSpacing(12);
        tab->setLayout(layout);

        auto* titleLabel = new QLabel(title);
        QFont tf = titleLabel->font();
        tf.setBold(true);
        tf.setPointSize(std::max(11, tf.pointSize() + 2));
        titleLabel->setFont(tf);
        layout->addWidget(titleLabel);

        auto* hint = new QLabel(desc);
        hint->setProperty("hint", true);
        hint->setWordWrap(true);
        layout->addWidget(hint);

        auto* box = new QGroupBox("Planned");
        auto* boxLayout = new QVBoxLayout;
        auto* list = new QLabel;
        list->setTextFormat(Qt::RichText);
        list->setWordWrap(true);
        QString html = "<ul>";
        for (const auto& item : bullets) html += "<li>" + item.toHtmlEscaped() + "</li>";
        html += "</ul>";
        list->setText(html);
        boxLayout->addWidget(list);
        box->setLayout(boxLayout);

        layout->addWidget(box);
        layout->addStretch(1);
    };

    setupRoadmapTab(
        manufacturingTab,
        "Manufacturing (Additive & Coating)",
        "Planned presets and helpers for process-driven experiments. This page is a placeholder.",
        QStringList()
            << "Additive manufacturing presets (scan path / layer schedule)"
            << "Coating / deposition parameter templates"
            << "Import simple toolpath formats (CSV / G-code) for boundary conditions"
            << "Run + compare multiple process recipes");

    setupRoadmapTab(
        transformationTab,
        "Structural Transformation",
        "Planned analysis for phase/structure evolution across checkpoints. This page is a placeholder.",
        QStringList()
            << "Track order-parameter statistics over time"
            << "Compare checkpoint snapshots (Δ images / metrics)"
            << "Optional peak-based orientation/defect analysis");

    setupRoadmapTab(
        mechanicsTab,
        "Stress–Strain",
        "Planned mechanical loading helpers. This page is a placeholder.",
        QStringList()
            << "Uniaxial / shear loading presets (small-scale demos)"
            << "Stress–strain curve export (CSV + plots)"
            << "Batch sweep on loading rate / amplitude");

    refreshResultsFileList();

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

    connectValueChanged(sweepU0Single_);
    connectValueChanged(sweepU0Start_);
    connectValueChanged(sweepU0End_);
    connectValueChanged(sweepU0Step_);
    connectValueChanged(sweepCon0Single_);
    connectValueChanged(sweepCon0Start_);
    connectValueChanged(sweepCon0End_);
    connectValueChanged(sweepCon0Step_);
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

void MainWindow::setResultsDir(const QString& dirPath) {
    if (!resultsDir_) return;
    resultsDir_->setText(dirPath);
    refreshResultsFileList();
}

void MainWindow::refreshResultsFileList() {
    if (!resultsDir_ || !resultsFiles_ || !resultsInfo_ || !resultsPlot_ || !resultsYColumn_ || !resultsTable_) return;

    const QString prev = resultsFiles_->currentItem() ? resultsFiles_->currentItem()->text() : QString();

    resultsFiles_->clear();
    resultsCurrentFile_.clear();
    resultsImageEnabled_ = false;
    resultsColumns_.clear();
    resultsPlot_->clear();
    resultsInfo_->setText("Pick a file to visualize.");
    resultsYColumn_->clear();
    resultsYColumn_->setEnabled(false);
    resultsTable_->clear();
    resultsTable_->setRowCount(0);
    resultsTable_->setColumnCount(0);
    if (resultsImageValue_) {
        resultsImageValue_->clear();
        resultsImageValue_->setEnabled(false);
    }
    if (resultsImagePointSize_) {
        resultsImagePointSize_->setEnabled(false);
    }
    if (resultsImageInfo_) {
        resultsImageInfo_->setText("Select a Phimax_*.txt snapshot file to render.");
    }
    if (resultsImageItem_) {
        resultsImageItem_->setPixmap(QPixmap());
    }
    if (resultsImageView_) {
        resultsImageView_->resetTransform();
    }

    const QString dirPath = resultsDir_->text().trimmed();
    if (dirPath.isEmpty()) {
        resultsInfo_->setText("Run directory is empty.");
        return;
    }
    QDir dir(dirPath);
    if (!dir.exists()) {
        resultsInfo_->setText("Directory not found: " + dirPath);
        return;
    }

    const QStringList files = dir.entryList(QStringList() << "*.txt", QDir::Files, QDir::Name);
    if (files.isEmpty()) {
        resultsInfo_->setText("No .txt files found in: " + dirPath);
        return;
    }

    for (const QString& file : files) {
        auto* item = new QListWidgetItem(file);
        item->setToolTip(dir.filePath(file));
        resultsFiles_->addItem(item);
    }

    resultsInfo_->setText(QString("Found %1 .txt file(s). Select one to visualize.").arg(files.size()));

    if (!prev.isEmpty()) {
        const QList<QListWidgetItem*> matches = resultsFiles_->findItems(prev, Qt::MatchExactly);
        if (!matches.isEmpty()) {
            resultsFiles_->setCurrentItem(matches.first());
            return;
        }
    }
    resultsFiles_->setCurrentRow(0);
}

static QVector<double> extractNumbersFromLine(const QString& line) {
    static const QRegularExpression re(R"([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)");
    QVector<double> nums;
    auto it = re.globalMatch(line);
    while (it.hasNext()) {
        const auto m = it.next();
        bool ok = false;
        const double v = m.captured(0).toDouble(&ok);
        if (ok) nums.push_back(v);
    }
    return nums;
}

void MainWindow::loadResultsFile(const QString& filePath) {
    if (!resultsInfo_ || !resultsPlot_ || !resultsYColumn_ || !resultsTable_) return;

    resultsCurrentFile_ = filePath;
    resultsImageEnabled_ = false;
    resultsColumns_.clear();
    resultsPlot_->clear();
    resultsYColumn_->clear();
    resultsYColumn_->setEnabled(false);
    resultsTable_->clear();
    resultsTable_->setRowCount(0);
    resultsTable_->setColumnCount(0);
    if (resultsImageValue_) {
        resultsImageValue_->clear();
        resultsImageValue_->setEnabled(false);
    }
    if (resultsImagePointSize_) {
        resultsImagePointSize_->setEnabled(false);
    }
    if (resultsImageInfo_) {
        resultsImageInfo_->setText("Select a Phimax_*.txt snapshot file to render.");
    }
    if (resultsImageItem_) {
        resultsImageItem_->setPixmap(QPixmap());
    }

    QFile f(filePath);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        resultsInfo_->setText("Failed to open: " + filePath);
        return;
    }

    const QString name = QFileInfo(filePath).fileName();
    const QString base = name.toLower();
    const bool isPhimaxSnapshot = base.startsWith("phimax_") || base.startsWith("phimaxq4q6_") || base.startsWith("phimaxq");

    bool skipKommentarHeader = false;
    {
        QTextStream probe(&f);
        QString first;
        QString second;
        while (!probe.atEnd()) {
            const QString line = probe.readLine().trimmed();
            if (line.isEmpty() || line.startsWith('#')) continue;
            if (first.isEmpty()) {
                first = line;
                continue;
            }
            second = line;
            break;
        }
        static const QRegularExpression intOnly(R"(^[+-]?\d+$)");
        skipKommentarHeader = intOnly.match(first).hasMatch() && second.toLower().contains("kommentar");
        f.seek(0);
    }

    QTextStream in(&f);
    QVector<QVector<double>> rows;
    rows.reserve(2048);
    int maxCols = 0;
    int meaningful = 0;
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;
        if (line.startsWith('#')) continue;
        meaningful++;
        if (skipKommentarHeader && meaningful <= 2) continue;
        const QVector<double> nums = extractNumbersFromLine(line);
        if (nums.isEmpty()) continue;
        maxCols = std::max(maxCols, static_cast<int>(nums.size()));
        rows.push_back(nums);
    }

    if (rows.isEmpty() || maxCols <= 0) {
        resultsInfo_->setText("No numeric data found in: " + name);
        return;
    }

    resultsColumns_.resize(maxCols);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto& row : rows) {
        for (int c = 0; c < maxCols; ++c) {
            resultsColumns_[c].push_back(c < row.size() ? row[c] : nan);
        }
    }

    resultsTable_->setRowCount(rows.size());
    resultsTable_->setColumnCount(maxCols);
    QStringList headers;
    headers.reserve(maxCols);
    for (int c = 0; c < maxCols; ++c) headers << QString("col%1").arg(c);
    resultsTable_->setHorizontalHeaderLabels(headers);
    for (int r = 0; r < rows.size(); ++r) {
        for (int c = 0; c < maxCols; ++c) {
            const double v = (c < rows[r].size()) ? rows[r][c] : nan;
            auto* item = new QTableWidgetItem(std::isfinite(v) ? QString::number(v, 'g', 12) : QString());
            resultsTable_->setItem(r, c, item);
        }
    }

    auto colLabel = [&](int c) -> QString {
        if (maxCols <= 1) return "value";
        if (base.contains("energy")) {
            if (c == 0) return "step";
            if (c == 1) return "energy";
            if (c == 2) return "density";
        }
        if (base == "phimax.txt") {
            if (c == 0) return "step";
            if (c == 1) return "phic";
            if (c == 2) return "phimax";
            if (c == 3) return "myid";
        }
        if (isPhimaxSnapshot) {
            if (c == 0) return "x";
            if (c == 1) return "y";
            if (c == 2) return "z";
            if (c == 3) return "phimax";
            if (c == 4) return "con";
            if (c == 5) return "q4";
            if (c == 6) return "q6";
        }
        if (base.contains("checkpoint")) {
            if (c == 0) return "step";
            if (c == 1) return "timestamp_ms";
        }
        return QString("col%1").arg(c);
    };

    if (isPhimaxSnapshot) {
        if (resultsImageValue_) {
            resultsImageValue_->clear();
            if (maxCols >= 4) resultsImageValue_->addItem("phimax (col3)", 3);
            if (maxCols >= 5) resultsImageValue_->addItem("con (col4)", 4);
            if (maxCols >= 6) resultsImageValue_->addItem("q4 (col5)", 5);
            if (maxCols >= 7) resultsImageValue_->addItem("q6 (col6)", 6);
            resultsImageValue_->setEnabled(resultsImageValue_->count() > 0);
        }
        if (resultsImagePointSize_) {
            resultsImagePointSize_->setEnabled(true);
        }
        resultsImageEnabled_ = resultsImageValue_ && (resultsImageValue_->count() > 0);
        resultsInfo_->setText(QString("%1 • %2 row(s) • %3 col(s) • use Image tab").arg(name).arg(rows.size()).arg(maxCols));
        if (resultsViewTabs_ && resultsImagePage_) resultsViewTabs_->setCurrentWidget(resultsImagePage_);
        renderResultsImage();
        return;
    }

    resultsYColumn_->clear();
    if (maxCols <= 1) {
        resultsYColumn_->addItem(colLabel(0), 0);
        resultsYColumn_->setEnabled(false);
        resultsInfo_->setText(QString("%1 • %2 row(s) • %3 col(s)").arg(name).arg(rows.size()).arg(maxCols));
        applyResultsYColumn(0);
        return;
    }

    for (int c = 1; c < maxCols; ++c) {
        resultsYColumn_->addItem(colLabel(c) + QString(" (col%1)").arg(c), c);
    }
    resultsYColumn_->setEnabled(true);
    resultsYColumn_->setCurrentIndex(0);

    resultsInfo_->setText(QString("%1 • %2 row(s) • %3 col(s) • x=%4").arg(name).arg(rows.size()).arg(maxCols).arg(colLabel(0)));
    applyResultsYColumn(resultsYColumn_->currentData().toInt());
}

void MainWindow::applyResultsYColumn(int columnIndex) {
    if (!resultsPlot_ || resultsColumns_.isEmpty()) return;

    const int cols = resultsColumns_.size();
    const int n = resultsColumns_.at(0).size();
    if (n <= 0) {
        resultsPlot_->clear();
        return;
    }

    QVector<double> xs;
    QVector<double> ys;
    xs.reserve(n);
    ys.reserve(n);

    QString yLabel;
    if (cols <= 1) {
        for (int i = 0; i < n; ++i) {
            xs.push_back(i);
            ys.push_back(resultsColumns_.at(0).at(i));
        }
        yLabel = "value";
    } else {
        const int ycol = std::clamp(columnIndex, 1, cols - 1);
        for (int i = 0; i < n; ++i) {
            xs.push_back(resultsColumns_.at(0).at(i));
            ys.push_back(resultsColumns_.at(ycol).at(i));
        }
        yLabel = QString("col%1").arg(ycol);
    }

    resultsPlot_->setData(std::move(xs), std::move(ys), yLabel);
}

static QColor jetLikeColor(double t) {
    t = std::clamp(t, 0.0, 1.0);
    // A compact "jet-like" ramp that matches common phase-field screenshots:
    // deep blue -> cyan/green -> yellow -> bright red.
    const double u = t * 0.875;  // keep the red endpoint bright
    auto ch = [](double x) {
        return std::clamp(1.5 - std::abs(x), 0.0, 1.0);
    };
    const double r = ch(4.0 * u - 3.0);
    const double g = ch(4.0 * u - 2.0);
    const double b = ch(4.0 * u - 1.0);
    return QColor::fromRgbF(r, g, b);
}

// ============================================================================
// Heatmap helper functions
// ============================================================================

/// @brief Pre-computed 256-level color lookup table for jet-like colormap.
/// @return Static array of QRgb colors (blue -> red gradient).
static const std::array<QRgb, 256>& jetLikeLut() {
    static const std::array<QRgb, 256> lut = []() {
        std::array<QRgb, 256> out{};
        for (int i = 0; i < static_cast<int>(out.size()); ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(out.size() - 1);
            out[static_cast<size_t>(i)] = jetLikeColor(t).rgba();
        }
        return out;
    }();
    return lut;
}

/// @brief Smooth Hermite interpolation (3rd order).
/// @param edge0 Lower edge of the transition range.
/// @param edge1 Upper edge of the transition range.
/// @param x Input value to interpolate.
/// @return Smoothstep result in [0, 1].
static double smoothstep(double edge0, double edge1, double x) {
    if (edge0 == edge1) return x < edge0 ? 0.0 : 1.0;
    const double t = std::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

/// @brief Linearly interpolate two QRgb colors.
/// @param a Start color (t=0).
/// @param b End color (t=1).
/// @param t Interpolation factor in [0, 1].
/// @return Interpolated QRgb color.
static QRgb lerpRgb(QRgb a, QRgb b, double t) {
    t = std::clamp(t, 0.0, 1.0);
    const int ar = qRed(a), ag = qGreen(a), ab = qBlue(a);
    const int br = qRed(b), bg = qGreen(b), bb = qBlue(b);
    const int r = std::clamp(static_cast<int>(std::lround(ar + (br - ar) * t)), 0, 255);
    const int g = std::clamp(static_cast<int>(std::lround(ag + (bg - ag) * t)), 0, 255);
    const int bl = std::clamp(static_cast<int>(std::lround(ab + (bb - ab) * t)), 0, 255);
    return qRgb(r, g, bl);
}

/// @brief Extract bit pattern of a double (for precise hashing).
/// @param v Double value.
/// @return 64-bit representation, with -0.0 normalized to +0.0.
static quint64 doubleBits(double v) {
    quint64 bits = 0;
    static_assert(sizeof(bits) == sizeof(v), "unexpected double size");
    std::memcpy(&bits, &v, sizeof(bits));
    // Normalize -0.0 to +0.0 for deduplication.
    if (bits == 0x8000000000000000ULL) bits = 0;
    return bits;
}

/// @brief Key for hashing (x, y) coordinate pairs by bit pattern.
struct XYKey {
    quint64 xb = 0;
    quint64 yb = 0;
};

static inline bool operator==(const XYKey& a, const XYKey& b) noexcept {
    return a.xb == b.xb && a.yb == b.yb;
}

/// @brief Hash function for XYKey (required by QHash).
inline uint qHash(const XYKey& key, uint seed = 0) noexcept {
    // Boost-like 64-bit hash combining.
    const quint64 mixed = key.xb ^ (key.yb + 0x9e3779b97f4a7c15ULL + (key.xb << 6) + (key.xb >> 2));
    return qHash(mixed, seed);
}

void MainWindow::renderResultsImage() {
    if (!resultsImageEnabled_ || !resultsImageItem_ || !resultsImageScene_ || !resultsImageView_ || !resultsImageInfo_) {
        return;
    }
    if (!resultsImageValue_ || !resultsImagePointSize_) {
        resultsImageInfo_->setText("Image controls are not available.");
        return;
    }
    if (resultsImageValue_->count() <= 0) {
        resultsImageInfo_->setText("No value columns available for this file.");
        return;
    }
    if (resultsColumns_.size() < 3 || resultsColumns_.at(0).isEmpty()) {
        resultsImageInfo_->setText("No point data loaded.");
        return;
    }

    const int n = resultsColumns_.at(0).size();
    const int xCol = 0;
    const int yCol = 1;
    const int valueCol = resultsImageValue_->currentData().toInt();
    if (valueCol < 0 || valueCol >= resultsColumns_.size()) {
        resultsImageInfo_->setText("Invalid value column.");
        return;
    }

    // ========================================================================
    // Step 1: Deduplicate by exact (x, y) coordinates and average values
    // ========================================================================
    struct Accum {
        double x = 0.0;
        double y = 0.0;
        double sumV = 0.0;
        int count = 0;
    };
    QHash<XYKey, Accum> dedup;
    dedup.reserve(n);
    int rawCount = 0;

    for (int i = 0; i < n; ++i) {
        const double x = resultsColumns_.at(xCol).at(i);
        const double y = resultsColumns_.at(yCol).at(i);
        const double v = resultsColumns_.at(valueCol).at(i);
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(v)) continue;
        rawCount++;
        const XYKey key{doubleBits(x), doubleBits(y)};
        auto it = dedup.find(key);
        if (it == dedup.end()) {
            dedup.insert(key, Accum{x, y, v, 1});
        } else {
            it->sumV += v;
            it->count += 1;
        }
    }

    if (dedup.isEmpty()) {
        resultsImageInfo_->setText("No finite point data found.");
        resultsImageItem_->setPixmap(QPixmap());
        return;
    }

    // Build deduplicated point list with averaged values
    struct Pt {
        double x = 0.0;
        double y = 0.0;
        double v = 0.0;
    };
    QVector<Pt> points;
    points.reserve(dedup.size());

    double minX = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    double minV = std::numeric_limits<double>::infinity();
    double maxV = -std::numeric_limits<double>::infinity();

    for (auto it = dedup.constBegin(); it != dedup.constEnd(); ++it) {
        const Accum& a = it.value();
        if (a.count <= 0) continue;
        const double vAvg = a.sumV / static_cast<double>(a.count);
        points.push_back(Pt{a.x, a.y, vAvg});
        minX = std::min(minX, a.x);
        maxX = std::max(maxX, a.x);
        minY = std::min(minY, a.y);
        maxY = std::max(maxY, a.y);
        minV = std::min(minV, vAvg);
        maxV = std::max(maxV, vAvg);
    }

    const int uniqueCount = points.size();
    if (uniqueCount <= 0 || !std::isfinite(minX) || !std::isfinite(maxX) || !std::isfinite(minY) || !std::isfinite(maxY)) {
        resultsImageInfo_->setText("No finite point data found.");
        resultsImageItem_->setPixmap(QPixmap());
        return;
    }

    const double rangeX = std::max(1e-12, maxX - minX);
    const double rangeY = std::max(1e-12, maxY - minY);
    const int longSide = 1024;
    int width = longSide;
    int height = longSide;
    if (rangeX >= rangeY) {
        height = std::max(1, static_cast<int>(std::lround(static_cast<double>(longSide) * (rangeY / rangeX))));
    } else {
        width = std::max(1, static_cast<int>(std::lround(static_cast<double>(longSide) * (rangeX / rangeY))));
    }

    // ========================================================================
    // Step 2: Gaussian kernel interpolation
    // ========================================================================
    // Interpret "Smooth σ (px)" spinbox value as σ (Gaussian standard deviation).
    // Use search radius R = 3σ for Gaussian kernel: w(d) = exp(-d²/(2σ²))
    const int sigmaPx = std::clamp(resultsImagePointSize_->value(), 1, 64);
    const int radiusPx = sigmaPx * 3;  // Search radius R = 3σ
    const int r2 = radiusPx * radiusPx;
    const double sigma = std::max(1e-6, static_cast<double>(sigmaPx));
    const double invTwoSigma2 = 1.0 / (2.0 * sigma * sigma);

    // Pre-compute Gaussian weights by squared distance
    QVector<float> wByD2(r2 + 1, 0.0f);
    for (int d2 = 0; d2 <= r2; ++d2) {
        wByD2[d2] = static_cast<float>(std::exp(-static_cast<double>(d2) * invTwoSigma2));
    }

    // Reuse buffers if possible
    const int pixels = width * height;
    const QSize newSize(width, height);
    if (resultsImageBufferSize_ != newSize) {
        resultsImageBufferSize_ = newSize;
        resultsImageSum_.resize(pixels);
        resultsImageWsum_.resize(pixels);
        resultsImageWmax_.resize(pixels);
    }
    // Clear buffers
    std::fill(resultsImageSum_.begin(), resultsImageSum_.end(), 0.0f);
    std::fill(resultsImageWsum_.begin(), resultsImageWsum_.end(), 0.0f);
    std::fill(resultsImageWmax_.begin(), resultsImageWmax_.end(), 0.0f);

    const double sx = (width - 1) / rangeX;
    const double sy = (height - 1) / rangeY;

    // Splat each point into its neighborhood
    for (const Pt& pt : points) {
        const int cx = std::clamp(static_cast<int>(std::lround((pt.x - minX) * sx)), 0, width - 1);
        const int cy = std::clamp(static_cast<int>(std::lround((maxY - pt.y) * sy)), 0, height - 1);
        const int x0 = std::max(0, cx - radiusPx);
        const int x1 = std::min(width - 1, cx + radiusPx);
        const int y0 = std::max(0, cy - radiusPx);
        const int y1 = std::min(height - 1, cy + radiusPx);

        const float v = static_cast<float>(pt.v);
        for (int py = y0; py <= y1; ++py) {
            const int dy = py - cy;
            const int dy2 = dy * dy;
            const int rowOffset = py * width;
            for (int px = x0; px <= x1; ++px) {
                const int dx = px - cx;
                const int d2 = dx * dx + dy2;
                if (d2 > r2) continue;
                const float w = wByD2[d2];
                const int idx = rowOffset + px;
                resultsImageWsum_[idx] += w;
                resultsImageSum_[idx] += w * v;
                resultsImageWmax_[idx] = std::max(resultsImageWmax_[idx], w);
            }
        }
    }

    // ========================================================================
    // Step 3: Colorize with boundary whitespace
    // ========================================================================
    // wMax ≈ exp(-d_min²/(2σ²)) gives coverage mask.
    // Start fading to white at ~3σ, fully colored by ~2.5σ.
    constexpr double kFadeStartSigma = 3.0;
    constexpr double kFadeEndSigma = 2.5;
    const double fadeLo = std::exp(-0.5 * kFadeStartSigma * kFadeStartSigma);
    const double fadeHi = std::exp(-0.5 * kFadeEndSigma * kFadeEndSigma);

    const auto& lut = jetLikeLut();
    const QRgb white = qRgb(255, 255, 255);
    const double invRangeV = (std::isfinite(minV) && std::isfinite(maxV) && (maxV - minV) > 1e-12)
                                 ? (1.0 / (maxV - minV))
                                 : 0.0;

    QImage img(width, height, QImage::Format_ARGB32);
    img.fill(Qt::white);

    for (int y = 0; y < height; ++y) {
        auto* dst = reinterpret_cast<QRgb*>(img.scanLine(y));
        const int row = y * width;
        for (int x = 0; x < width; ++x) {
            const int idx = row + x;
            const float wSum = resultsImageWsum_[idx];
            if (wSum <= 0.0f) {
                dst[x] = white;
                continue;
            }
            const float wMax = resultsImageWmax_[idx];
            if (wMax <= static_cast<float>(fadeLo)) {
                dst[x] = white;
                continue;
            }

            // Normalized field value
            const double v = static_cast<double>(resultsImageSum_[idx]) / static_cast<double>(wSum);
            double t = invRangeV > 0.0 ? ((v - minV) * invRangeV) : 0.5;
            t = std::clamp(t, 0.0, 1.0);
            const int li = std::clamp(static_cast<int>(std::lround(t * 255.0)), 0, 255);
            const QRgb c = lut[static_cast<size_t>(li)];

            // Fade to white if coverage is low
            const double alpha = smoothstep(fadeLo, fadeHi, static_cast<double>(wMax));
            dst[x] = lerpRgb(white, c, alpha);
        }
    }

    resultsImageItem_->setPixmap(QPixmap::fromImage(img));
    resultsImageScene_->setSceneRect(QRectF(QPointF(0, 0), QSizeF(width, height)));

    const QString valueLabel = resultsImageValue_->currentText();
    resultsImageInfo_->setText(
        QString("%1 • points=%2 • unique(x,y)=%3 • σ=%4px • x=[%5,%6] y=[%7,%8] • %9=[%10,%11] • Ctrl+wheel zoom")
            .arg(QFileInfo(resultsCurrentFile_).fileName())
            .arg(rawCount)
            .arg(uniqueCount)
            .arg(sigmaPx)
            .arg(minX, 0, 'g', 4)
            .arg(maxX, 0, 'g', 4)
            .arg(minY, 0, 'g', 4)
            .arg(maxY, 0, 'g', 4)
            .arg(valueLabel)
            .arg(minV, 0, 'g', 4)
            .arg(maxV, 0, 'g', 4));

    fitResultsImage();
}

void MainWindow::fitResultsImage() {
    if (!resultsImageView_ || !resultsImageItem_) return;
    if (resultsImageItem_->pixmap().isNull()) return;
    resultsImageView_->resetTransform();
    resultsImageView_->fitInView(resultsImageItem_, Qt::KeepAspectRatio);
}

void MainWindow::exportResultsImage() {
    if (!resultsImageItem_ || resultsImageItem_->pixmap().isNull()) {
        QMessageBox::information(this, "Export image", "No image to export. Select a Phimax_*.txt file first.");
        return;
    }

    const QString baseDir = resultsDir_ ? resultsDir_->text().trimmed() : QString();
    const QString suggestDir = baseDir.isEmpty() ? QStandardPaths::writableLocation(QStandardPaths::PicturesLocation) : baseDir;
    const QString baseName = QFileInfo(resultsCurrentFile_).completeBaseName();
    const QString suggested = QDir(suggestDir).filePath(baseName + ".png");

    const QString fileName = QFileDialog::getSaveFileName(this, "Export PNG", suggested, "PNG Image (*.png)");
    if (fileName.isEmpty()) return;

    if (!resultsImageItem_->pixmap().toImage().save(fileName)) {
        QMessageBox::warning(this, "Export image", "Failed to save: " + fileName);
    }
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

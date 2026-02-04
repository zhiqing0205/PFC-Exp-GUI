#include "MainWindow.h"

#include <QCheckBox>
#include <QComboBox>
#include <QCoreApplication>
#include <QDateTime>
#include <QDesktopServices>
#include <QDialog>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QFrame>
#include <QDialogButtonBox>
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
#include <QMouseEvent>
#include <QPainter>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QRegularExpression>
#include <QTextStream>
#include <QScrollArea>
#include <QScrollBar>
#include <QSplitter>
#include <QSpinBox>
#include <QStackedWidget>
#include <QStandardPaths>
#include <QStatusBar>
#include <QStyle>
#include <QSysInfo>
#include <QTabBar>
#include <QTabWidget>
#include <QSignalBlocker>
#include <QToolButton>
#include <QUrl>
#include <QVBoxLayout>
#include <QFontDatabase>
#include <QTableWidget>
#include <QHeaderView>
#include <QWheelEvent>
#include <QResizeEvent>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>

class NoWheelSpinBox final : public QSpinBox {
public:
    using QSpinBox::QSpinBox;

protected:
    void wheelEvent(QWheelEvent* event) override {
        // Prevent accidental value changes while scrolling the page.
        event->ignore();
    }
};

class NoWheelDoubleSpinBox final : public QDoubleSpinBox {
public:
    using QDoubleSpinBox::QDoubleSpinBox;

protected:
    void wheelEvent(QWheelEvent* event) override {
        // Prevent accidental value changes while scrolling the page.
        event->ignore();
    }
};

static QDoubleSpinBox* makeDoubleSpin(double min, double max, double value, int decimals, double step, bool disableWheel = false) {
    auto* box = disableWheel ? static_cast<QDoubleSpinBox*>(new NoWheelDoubleSpinBox) : new QDoubleSpinBox;
    box->setRange(min, max);
    box->setDecimals(decimals);
    box->setSingleStep(step);
    box->setValue(value);
    box->setKeyboardTracking(false);
    return box;
}

static QSpinBox* makeIntSpin(int min, int max, int value, int step, bool disableWheel = false) {
    auto* box = disableWheel ? static_cast<QSpinBox*>(new NoWheelSpinBox) : new QSpinBox;
    box->setRange(min, max);
    box->setSingleStep(step);
    box->setValue(value);
    box->setKeyboardTracking(false);
    return box;
}

static QString pickExistingDirectory(QWidget* parent, const QString& caption, const QString& baseDir) {
    QFileDialog dlg(parent, caption, baseDir);
    dlg.setFileMode(QFileDialog::Directory);
    dlg.setOption(QFileDialog::ShowDirsOnly, true);
    dlg.setOption(QFileDialog::DontUseNativeDialog, true);
    if (dlg.exec() != QDialog::Accepted) return QString();
    const QStringList selected = dlg.selectedFiles();
    if (selected.isEmpty()) return QString();
    return selected.first();
}

static QString pickSavePngFile(QWidget* parent, const QString& caption, const QString& suggestedPath) {
    QFileDialog dlg(parent, caption, suggestedPath);
    dlg.setAcceptMode(QFileDialog::AcceptSave);
    dlg.setNameFilter("PNG Image (*.png)");
    dlg.setDefaultSuffix("png");
    dlg.setOption(QFileDialog::DontUseNativeDialog, true);
    if (dlg.exec() != QDialog::Accepted) return QString();
    const QStringList selected = dlg.selectedFiles();
    if (selected.isEmpty()) return QString();
    return selected.first();
}

class ZoomableGraphicsView final : public QGraphicsView {
public:
    using QGraphicsView::QGraphicsView;

protected:
    void wheelEvent(QWheelEvent* event) override {
        int delta = event->angleDelta().y();
        if (delta == 0) delta = event->pixelDelta().y();
        if (delta == 0) {
            QGraphicsView::wheelEvent(event);
            return;
        }

        constexpr double base = 1.0015;
        double factor = std::pow(base, static_cast<double>(delta));

        constexpr double kMinScale = 0.05;
        constexpr double kMaxScale = 200.0;
        const double current = transform().m11();
        const double next = current * factor;
        if (next < kMinScale) factor = kMinScale / current;
        if (next > kMaxScale) factor = kMaxScale / current;

        scale(factor, factor);
        event->accept();
    }

    void mousePressEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton && dragMode() == QGraphicsView::ScrollHandDrag) {
            setCursor(Qt::ClosedHandCursor);
        }
        QGraphicsView::mousePressEvent(event);
    }

    void mouseReleaseEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton && dragMode() == QGraphicsView::ScrollHandDrag) {
            setCursor(Qt::OpenHandCursor);
        }
        QGraphicsView::mouseReleaseEvent(event);
    }
};

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("MID Nano");
    resize(1100, 920);

    auto* uploadLicenseBtn = new QPushButton("Upload License…");
    uploadLicenseBtn->setObjectName("tabUploadLicense");
    uploadLicenseBtn->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    auto* aboutBtn = new QPushButton("About");
    aboutBtn->setObjectName("tabAbout");
    aboutBtn->setIcon(style()->standardIcon(QStyle::SP_MessageBoxInformation));

    auto* tabs = new QTabWidget;
    tabs->setObjectName("mainTabs");
    tabs->setDocumentMode(true);
    tabs->setUsesScrollButtons(true);
    tabs->tabBar()->setDrawBase(false);
    tabs->tabBar()->setExpanding(false);

    auto* tabActions = new QWidget;
    tabActions->setObjectName("mainTabsActions");
    tabActions->setAttribute(Qt::WA_StyledBackground, true);
    auto* tabActionsLayout = new QHBoxLayout;
    tabActionsLayout->setContentsMargins(8, 4, 8, 4);
    tabActionsLayout->setSpacing(6);
    tabActionsLayout->setAlignment(Qt::AlignVCenter);
    tabActionsLayout->addWidget(uploadLicenseBtn);
    tabActionsLayout->addWidget(aboutBtn);
    tabActions->setLayout(tabActionsLayout);
    tabs->setCornerWidget(tabActions, Qt::TopRightCorner);

    auto* expTab = new QWidget;
    auto* resultsTab = new QWidget;
    auto* manufacturingTab = new QWidget;
    auto* transformationTab = new QWidget;
    auto* mechanicsTab = new QWidget;

    auto* expScroll = new QScrollArea;
    expScroll->setFrameShape(QFrame::NoFrame);
    expScroll->setWidgetResizable(true);
    expScroll->setWidget(expTab);

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

    tabs->addTab(expScroll, "Experiment");
    tabs->addTab(resultsScroll, "Visualizer");
    tabs->addTab(manufacturingScroll, "Manufacturing");
    tabs->addTab(transformationScroll, "Transformation");
    tabs->addTab(mechanicsScroll, "Stress–Strain");

    connect(tabs, &QTabWidget::currentChanged, this, [this, tabs, resultsScroll](int) {
        if (!resultsScroll) return;
        if (tabs->currentWidget() != resultsScroll) return;
        if (!resultsDir_ || !resultsFiles_) return;

        const QString current = resultsDir_->text().trimmed();
        if (!currentOutputDir_.isEmpty() && current != currentOutputDir_) {
            setResultsDir(currentOutputDir_);
            return;
        }
        refreshResultsFileList();
    });

    log_ = new QPlainTextEdit;
    log_->setReadOnly(true);
    log_->setMaximumBlockCount(5000);
    {
        QFont f = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        if (f.pointSizeF() > 0) {
            f.setPointSizeF(std::max(10.0, f.pointSizeF()));
        } else {
            f.setPointSize(11);
        }
        log_->setFont(f);
    }

    logDialog_ = new QDialog(this);
    logDialog_->setWindowTitle("Log");
    logDialog_->resize(980, 620);
    {
        auto* dlgLayout = new QVBoxLayout;
        dlgLayout->setContentsMargins(12, 12, 12, 12);
        dlgLayout->setSpacing(10);
        dlgLayout->addWidget(log_);
        logDialog_->setLayout(dlgLayout);
    }

    auto* central = new QWidget;
    auto* centralLayout = new QVBoxLayout;
    centralLayout->setContentsMargins(14, 14, 14, 14);
    centralLayout->setSpacing(0);
    centralLayout->addWidget(tabs, 1);
    central->setLayout(centralLayout);
    setCentralWidget(central);

    auto* sb = new QStatusBar;
    sb->setSizeGripEnabled(true);
    setStatusBar(sb);

    logStatusButton_ = new QToolButton;
    logStatusButton_->setAutoRaise(true);
    logStatusButton_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    logStatusButton_->setIcon(style()->standardIcon(QStyle::SP_FileDialogDetailedView));
    logStatusButton_->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
    logStatusButton_->setMinimumWidth(0);
    statusLogText_ = "Ready (click to open log)";
    statusLogTooltip_ = statusLogText_;
    updateStatusLogText();
    sb->addWidget(logStatusButton_, 1);
    connect(logStatusButton_, &QToolButton::clicked, this, [this]() {
        if (!logDialog_) return;
        logDialog_->show();
        logDialog_->raise();
        logDialog_->activateWindow();
    });

    connect(aboutBtn, &QPushButton::clicked, this, &MainWindow::showAboutDialog);
    connect(uploadLicenseBtn, &QPushButton::clicked, this, &MainWindow::uploadLicense);

    // ---------- Experiment tab ----------
    auto* expLayout = new QVBoxLayout;
    expLayout->setSpacing(10);
    expTab->setLayout(expLayout);

    auto header = [](const QString& text) {
        auto* l = new QLabel(text);
        QFont f = l->font();
        f.setBold(true);
        l->setFont(f);
        l->setAlignment(Qt::AlignCenter);
        return l;
    };

    auto* expParamsBox = new QGroupBox("Parameters");
    auto* expParamsGrid = new QGridLayout;
    expParamsGrid->setVerticalSpacing(6);
    expParamsGrid->setHorizontalSpacing(12);
    expParamsGrid->setColumnStretch(2, 1);

    expParamsGrid->addWidget(header("Param"), 0, 0);
    expParamsGrid->addWidget(header("Mode"), 0, 1);
    expParamsGrid->addWidget(header("Values"), 0, 2);

    sweepParams_.clear();
    sweepParams_.reserve(16);

    auto connectSpinChanged = [this](QObject* w) {
        if (!w) return;
        if (auto* d = qobject_cast<QDoubleSpinBox*>(w)) {
            connect(d, &QDoubleSpinBox::valueChanged, this, &MainWindow::updateExperimentPreview);
        } else if (auto* s = qobject_cast<QSpinBox*>(w)) {
            connect(s, &QSpinBox::valueChanged, this, &MainWindow::updateExperimentPreview);
        } else if (auto* e = qobject_cast<QLineEdit*>(w)) {
            connect(e, &QLineEdit::textChanged, this, &MainWindow::updateExperimentPreview);
        } else if (auto* c = qobject_cast<QComboBox*>(w)) {
            connect(c, &QComboBox::currentIndexChanged, this, &MainWindow::updateExperimentPreview);
        }
    };

    auto addDoubleParam = [&](const QString& key, const QString& labelText, double minV, double maxV, double defV, int decimals, double stepV) {
        SweepParamWidgets p;
        p.key = key;
        p.label = labelText;
        p.isInt = false;
        p.decimals = decimals;

        auto* name = new QLabel(labelText);
        name->setAlignment(Qt::AlignRight | Qt::AlignVCenter);

        p.mode = new QComboBox;
        p.mode->setProperty("kind", "mode");
        p.mode->addItem("Fixed", 0);
        p.mode->addItem("Range", 1);
        p.mode->addItem("List", 2);
        p.mode->setCurrentIndex(0);

        p.stack = new QStackedWidget;

        // Fixed page
        auto* fixedPage = new QWidget;
        auto* fixedLayout = new QHBoxLayout;
        fixedLayout->setContentsMargins(0, 0, 0, 0);
        p.fixedDouble = makeDoubleSpin(minV, maxV, defV, decimals, stepV, true);
        fixedLayout->addWidget(p.fixedDouble, 1);
        fixedPage->setLayout(fixedLayout);

        // Range page
        auto* rangePage = new QWidget;
        auto* rangeLayout = new QHBoxLayout;
        rangeLayout->setContentsMargins(0, 0, 0, 0);
        rangeLayout->setSpacing(8);
        rangeLayout->addWidget(new QLabel("Start"));
        p.rangeStartDouble = makeDoubleSpin(minV, maxV, defV, decimals, stepV, true);
        rangeLayout->addWidget(p.rangeStartDouble);
        rangeLayout->addWidget(new QLabel("End"));
        p.rangeEndDouble = makeDoubleSpin(minV, maxV, defV, decimals, stepV, true);
        rangeLayout->addWidget(p.rangeEndDouble);
        rangeLayout->addWidget(new QLabel("Step"));
        const double maxStep = std::max(1e-12, maxV - minV);
        p.rangeStepDouble = makeDoubleSpin(1e-12, maxStep, stepV, decimals, stepV, true);
        rangeLayout->addWidget(p.rangeStepDouble);
        rangeLayout->addWidget(new QLabel("List"));
        p.rangePreview = new QLineEdit;
        p.rangePreview->setReadOnly(true);
        p.rangePreview->setPlaceholderText("e.g. 0.0,0.1,0.2");
        rangeLayout->addWidget(p.rangePreview, 1);
        rangePage->setLayout(rangeLayout);

        // List page
        auto* listPage = new QWidget;
        auto* listLayout = new QHBoxLayout;
        listLayout->setContentsMargins(0, 0, 0, 0);
        auto* listLabel = new QLabel("Values");
        listLayout->addWidget(listLabel);
        p.listEdit = new QLineEdit;
        p.listEdit->setPlaceholderText("e.g. 0.05,0.10,0.15");
        listLayout->addWidget(p.listEdit, 1);
        listPage->setLayout(listLayout);

        p.stack->addWidget(fixedPage);
        p.stack->addWidget(rangePage);
        p.stack->addWidget(listPage);
        p.stack->setCurrentIndex(p.mode->currentIndex());

        connect(p.mode, &QComboBox::currentIndexChanged, this, [stack = p.stack](int idx) { stack->setCurrentIndex(idx); });

        // Trigger preview updates
        connectSpinChanged(p.mode);
        connectSpinChanged(p.fixedDouble);
        connectSpinChanged(p.rangeStartDouble);
        connectSpinChanged(p.rangeEndDouble);
        connectSpinChanged(p.rangeStepDouble);
        connectSpinChanged(p.listEdit);

        const int row = expParamsGrid->rowCount();
        expParamsGrid->addWidget(name, row, 0);
        expParamsGrid->addWidget(p.mode, row, 1);
        expParamsGrid->addWidget(p.stack, row, 2);
        sweepParams_.push_back(p);
    };

    auto addIntParam = [&](const QString& key, const QString& labelText, int minV, int maxV, int defV, int stepV, const QString& listPlaceholder) {
        SweepParamWidgets p;
        p.key = key;
        p.label = labelText;
        p.isInt = true;

        auto* name = new QLabel(labelText);
        name->setAlignment(Qt::AlignRight | Qt::AlignVCenter);

        p.mode = new QComboBox;
        p.mode->setProperty("kind", "mode");
        p.mode->addItem("Fixed", 0);
        p.mode->addItem("Range", 1);
        p.mode->addItem("List", 2);
        p.mode->setCurrentIndex(0);

        p.stack = new QStackedWidget;

        auto* fixedPage = new QWidget;
        auto* fixedLayout = new QHBoxLayout;
        fixedLayout->setContentsMargins(0, 0, 0, 0);
        p.fixedInt = makeIntSpin(minV, maxV, defV, stepV, true);
        fixedLayout->addWidget(p.fixedInt, 1);
        fixedPage->setLayout(fixedLayout);

        auto* rangePage = new QWidget;
        auto* rangeLayout = new QHBoxLayout;
        rangeLayout->setContentsMargins(0, 0, 0, 0);
        rangeLayout->setSpacing(8);
        rangeLayout->addWidget(new QLabel("Start"));
        p.rangeStartInt = makeIntSpin(minV, maxV, defV, stepV, true);
        rangeLayout->addWidget(p.rangeStartInt);
        rangeLayout->addWidget(new QLabel("End"));
        p.rangeEndInt = makeIntSpin(minV, maxV, defV, stepV, true);
        rangeLayout->addWidget(p.rangeEndInt);
        rangeLayout->addWidget(new QLabel("Step"));
        p.rangeStepInt = makeIntSpin(1, std::max(1, maxV - minV), stepV, 1, true);
        rangeLayout->addWidget(p.rangeStepInt);
        rangeLayout->addWidget(new QLabel("List"));
        p.rangePreview = new QLineEdit;
        p.rangePreview->setReadOnly(true);
        p.rangePreview->setPlaceholderText(listPlaceholder);
        rangeLayout->addWidget(p.rangePreview, 1);
        rangePage->setLayout(rangeLayout);

        auto* listPage = new QWidget;
        auto* listLayout = new QHBoxLayout;
        listLayout->setContentsMargins(0, 0, 0, 0);
        listLayout->addWidget(new QLabel("Values"));
        p.listEdit = new QLineEdit;
        p.listEdit->setPlaceholderText(listPlaceholder);
        listLayout->addWidget(p.listEdit, 1);
        listPage->setLayout(listLayout);

        p.stack->addWidget(fixedPage);
        p.stack->addWidget(rangePage);
        p.stack->addWidget(listPage);
        p.stack->setCurrentIndex(p.mode->currentIndex());

        connect(p.mode, &QComboBox::currentIndexChanged, this, [stack = p.stack](int idx) { stack->setCurrentIndex(idx); });

        connectSpinChanged(p.mode);
        connectSpinChanged(p.fixedInt);
        connectSpinChanged(p.rangeStartInt);
        connectSpinChanged(p.rangeEndInt);
        connectSpinChanged(p.rangeStepInt);
        connectSpinChanged(p.listEdit);

        const int row = expParamsGrid->rowCount();
        expParamsGrid->addWidget(name, row, 0);
        expParamsGrid->addWidget(p.mode, row, 1);
        expParamsGrid->addWidget(p.stack, row, 2);
        sweepParams_.push_back(p);
    };

    // Order matters for the Cartesian product (outer -> inner).
    addDoubleParam("u0", "u0", -2.0, 2.0, 0.05, 6, 0.01);
    addDoubleParam("con0", "con0", 0.0, 1.0, 0.2, 6, 0.01);
    addDoubleParam("sig", "sig", 0.0, 2.0, 0.05, 6, 0.01);
    addDoubleParam("dt", "dt", 1e-6, 1.0, 0.05, 6, 0.01);
    addDoubleParam("dx", "dx", 1e-6, 10.0, 0.125, 6, 0.01);

    addIntParam("steps", "steps", 1, 100000000, 200, 10, "e.g. 100,200,300");
    addIntParam("mod", "mod", 1, 100000000, 25, 1, "e.g. 25,50,100");
    addIntParam("seed", "seed", 0, 2147483647, 20200604, 1, "e.g. 1,2,3");

    expParamsBox->setLayout(expParamsGrid);

    auto* expOutBox = new QGroupBox("Output");
    auto* expOutForm = new QFormLayout;
    expBaseOutDir_ = new QLineEdit(defaultBaseOutputDir());
    expName_ = new QLineEdit;
    expName_->setPlaceholderText("Optional (default: timestamp)");
    auto* expBrowse = new QToolButton;
    expBrowse->setText("Browse…");
    expBrowse->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    expBrowse->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    connect(expBrowse, &QToolButton::clicked, this, &MainWindow::browseExperimentOutDir);

    auto* expBaseRow = new QWidget;
    auto* expBaseRowLayout = new QHBoxLayout;
    expBaseRowLayout->setContentsMargins(0, 0, 0, 0);
    expBaseRowLayout->addWidget(expBaseOutDir_, 1);
    expBaseRowLayout->addWidget(expBrowse);
    expBaseRow->setLayout(expBaseRowLayout);

    expOutForm->addRow("Base output dir", expBaseRow);
    expOutForm->addRow("Experiment folder name", expName_);
    expOutBox->setLayout(expOutForm);

    expPreview_ = new QLabel;
    expPreview_->setProperty("hint", true);
    expPreview_->setWordWrap(true);
    expStepProgress_ = new QProgressBar;
    expStepProgress_->setRange(0, 1);
    expStepProgress_->setValue(0);
    expStepProgress_->setFormat("Step %v / %m");
    expStepProgress_->setTextVisible(true);
    expProgress_ = new QProgressBar;
    expProgress_->setRange(0, 100);
    expProgress_->setValue(0);

    auto* expButtons = new QWidget;
    auto* expButtonsLayout = new QHBoxLayout;
    expButtonsLayout->setContentsMargins(0, 0, 0, 0);
    auto* expRun = new QPushButton("Run");
    auto* expStop = new QPushButton("Stop");
    auto* expOpenOut = new QPushButton("Open Output");
    expRun->setProperty("primary", true);
    expStop->setProperty("danger", true);
    expRun->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    expStop->setIcon(style()->standardIcon(QStyle::SP_MediaStop));
    expOpenOut->setIcon(style()->standardIcon(QStyle::SP_DirOpenIcon));
    connect(expRun, &QPushButton::clicked, this, &MainWindow::startExperiment);
    connect(expStop, &QPushButton::clicked, this, &MainWindow::stopRun);
    connect(expOpenOut, &QPushButton::clicked, this, &MainWindow::openCurrentOutputDir);
    expButtonsLayout->addWidget(expRun);
    expButtonsLayout->addWidget(expStop);
    expButtonsLayout->addStretch(1);
    expButtonsLayout->addWidget(expOpenOut);
    expButtons->setLayout(expButtonsLayout);

    expLayout->addWidget(expParamsBox);
    expLayout->addWidget(expPreview_);
    expLayout->addWidget(expOutBox);
    expLayout->addWidget(expStepProgress_);
    expLayout->addWidget(expProgress_);
    expLayout->addWidget(expButtons);
    expLayout->addStretch(1);

    // ---------- Results tab ----------
    auto* resultsLayout = new QVBoxLayout;
    resultsLayout->setSpacing(12);
    resultsTab->setLayout(resultsLayout);

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

    resultsStatus_ = new QLabel("Using current output directory when available.");
    resultsStatus_->setProperty("hint", true);
    resultsStatus_->setWordWrap(true);
    resultsLayout->addWidget(resultsStatus_);

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
    if (resultsViewTabs_->tabBar()) resultsViewTabs_->tabBar()->setDrawBase(false);
    resultsViewTabs_->setStyleSheet("QTabWidget::pane { border: 0px; }");

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

    auto* plotPage = new QWidget;
    auto* plotLayout = new QVBoxLayout;
    plotLayout->setContentsMargins(0, 0, 0, 0);
    plotLayout->setSpacing(10);

    auto* plotTopRow = new QWidget;
    auto* plotTopLayout = new QHBoxLayout;
    plotTopLayout->setContentsMargins(0, 0, 0, 0);
    plotTopLayout->addWidget(new QLabel("Value"));
    resultsImageValue_ = new QComboBox;
    resultsImageValue_->setMinimumWidth(180);
    resultsImageValue_->setEnabled(false);
    plotTopLayout->addWidget(resultsImageValue_);
    plotTopLayout->addWidget(new QLabel("Smooth σ (px)"));
    resultsImagePointSize_ = makeIntSpin(1, 64, 16, 1);
    resultsImagePointSize_->setToolTip("Gaussian smoothing radius (σ). Higher values create smoother heatmaps.");
    resultsImagePointSize_->setEnabled(false);
    plotTopLayout->addWidget(resultsImagePointSize_);
    plotTopLayout->addStretch(1);

    auto* imageFit = new QToolButton;
    imageFit->setText("Fit");
    imageFit->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    imageFit->setIcon(style()->standardIcon(QStyle::SP_DesktopIcon));
    plotTopLayout->addWidget(imageFit);

    auto* imageExport = new QToolButton;
    imageExport->setText("Export PNG…");
    imageExport->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    imageExport->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
    plotTopLayout->addWidget(imageExport);
    plotTopRow->setLayout(plotTopLayout);

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
    resultsImageView_->setCursor(Qt::OpenHandCursor);
    resultsImageView_->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    resultsImageView_->setResizeAnchor(QGraphicsView::AnchorUnderMouse);
    resultsImageItem_ = resultsImageScene_->addPixmap(QPixmap());

    plotLayout->addWidget(plotTopRow);
    plotLayout->addWidget(resultsImageView_, 1);
    plotLayout->addWidget(resultsImageInfo_);
    plotPage->setLayout(plotLayout);
    resultsPlotPage_ = plotPage;

    resultsViewTabs_->addTab(plotPage, "Plot");
    resultsViewTabs_->addTab(tablePage, "Table");
    const int plotTab = resultsViewTabs_->indexOf(plotPage);
    resultsViewTabs_->setTabEnabled(plotTab, false);
    if (resultsViewTabs_->tabBar()) resultsViewTabs_->tabBar()->setTabVisible(plotTab, false);
    resultsViewTabs_->setCurrentWidget(tablePage);

    resultsSplit->addWidget(filesBox);
    resultsSplit->addWidget(resultsViewTabs_);
    resultsSplit->setStretchFactor(0, 1);
    resultsSplit->setStretchFactor(1, 3);
    resultsSplit->setSizes(QList<int>() << 260 << 780);

    resultsLayout->addWidget(resultsSplit, 1);

    connect(resultsBrowse, &QToolButton::clicked, this, [this]() {
        appendLog("=== Visualizer: Browse clicked ===");
        const QString base = resultsDir_ ? resultsDir_->text() : defaultBaseOutputDir();
        const QString dir = pickExistingDirectory(this, "Select run directory", base);
        if (dir.isEmpty()) {
            appendLog("=== Visualizer: Browse cancelled ===");
            return;
        }
        appendLog("=== Visualizer: Browse selected: " + dir + " ===");
        setResultsDir(dir);
    });
    connect(resultsRefresh, &QToolButton::clicked, this, [this]() {
        appendLog("=== Visualizer: Refresh clicked ===");
        refreshResultsFileList();
    });
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

    connect(resultsImageValue_, &QComboBox::currentIndexChanged, this, [this](int) { renderResultsImage(); });
    connect(resultsImagePointSize_, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int) { renderResultsImage(); });
    connect(imageFit, &QToolButton::clicked, this, [this]() {
        appendLog("=== Visualizer: Fit clicked ===");
        fitResultsImage();
    });
    connect(imageExport, &QToolButton::clicked, this, [this]() {
        appendLog("=== Visualizer: Export PNG clicked ===");
        exportResultsImage();
    });

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
    if (expName_) connect(expName_, &QLineEdit::textChanged, this, &MainWindow::updateExperimentPreview);
    if (expBaseOutDir_) connect(expBaseOutDir_, &QLineEdit::textChanged, this, &MainWindow::updateExperimentPreview);

    updateExperimentPreview();
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
        return QDir(docs).filePath("MID Nano/outputs");
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
    const QString prev = resultsDir_->text().trimmed();
    resultsDir_->setText(dirPath);
    if (dirPath.trimmed() != prev) {
        resetResultsViewState();
        if (resultsFiles_) resultsFiles_->clear();
    }
    refreshResultsFileList();
}

void MainWindow::resetResultsViewState() {
    if (!resultsStatus_ || !resultsViewTabs_ || !resultsTable_) return;

    resultsCurrentFile_.clear();
    resultsImageEnabled_ = false;
    resultsColumns_.clear();

    resultsTable_->clear();
    resultsTable_->setRowCount(0);
    resultsTable_->setColumnCount(0);

    if (resultsPlotPage_) {
        const int plotTab = resultsViewTabs_->indexOf(resultsPlotPage_);
        if (plotTab >= 0) {
            resultsViewTabs_->setTabEnabled(plotTab, false);
            if (resultsViewTabs_->tabBar()) resultsViewTabs_->tabBar()->setTabVisible(plotTab, false);
        }
    }
    if (resultsViewTabs_->count() > 0) {
        const int plotTab = resultsPlotPage_ ? resultsViewTabs_->indexOf(resultsPlotPage_) : -1;
        const int tableTab = (plotTab == 0) ? 1 : 0;
        if (tableTab >= 0 && tableTab < resultsViewTabs_->count()) resultsViewTabs_->setCurrentIndex(tableTab);
    }
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
}

void MainWindow::refreshResultsFileList() {
    if (!resultsDir_ || !resultsStatus_ || !resultsFiles_ || !resultsViewTabs_ || !resultsTable_) return;

    const QString prev = resultsFiles_->currentItem() ? resultsFiles_->currentItem()->text() : QString();

    resultsFiles_->clear();
    resetResultsViewState();

    const QString dirPath = resultsDir_->text().trimmed();
    if (dirPath.isEmpty()) {
        resultsStatus_->setText("Run directory is empty.");
        return;
    }
    QDir dir(dirPath);
    if (!dir.exists()) {
        resultsStatus_->setText("Directory not found: " + dirPath);
        return;
    }

    const QStringList files = dir.entryList(QStringList() << "*.txt", QDir::Files, QDir::Name);
    if (files.isEmpty()) {
        resultsStatus_->setText("No .txt files found in: " + dirPath);
        return;
    }

    for (const QString& file : files) {
        auto* item = new QListWidgetItem(file);
        item->setToolTip(dir.filePath(file));
        resultsFiles_->addItem(item);
    }

    resultsStatus_->setText(QString("Found %1 .txt file(s). Select one to inspect.").arg(files.size()));

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
    if (!resultsStatus_ || !resultsViewTabs_ || !resultsTable_) return;

    resultsCurrentFile_ = filePath;
    resultsImageEnabled_ = false;
    resultsColumns_.clear();
    resultsTable_->clear();
    resultsTable_->setRowCount(0);
    resultsTable_->setColumnCount(0);
    if (auto* header = resultsTable_->horizontalHeader()) {
        header->setSectionResizeMode(QHeaderView::Interactive);
        header->setStretchLastSection(true);
    }
    if (resultsPlotPage_) {
        const int plotTab = resultsViewTabs_->indexOf(resultsPlotPage_);
        if (plotTab >= 0) {
            resultsViewTabs_->setTabEnabled(plotTab, false);
            if (resultsViewTabs_->tabBar()) resultsViewTabs_->tabBar()->setTabVisible(plotTab, false);
        }
    }
    if (resultsViewTabs_->count() > 0) {
        const int plotTab = resultsPlotPage_ ? resultsViewTabs_->indexOf(resultsPlotPage_) : -1;
        const int tableTab = (plotTab == 0) ? 1 : 0;
        if (tableTab >= 0 && tableTab < resultsViewTabs_->count()) resultsViewTabs_->setCurrentIndex(tableTab);
    }
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
        resultsStatus_->setText("Failed to open: " + filePath);
        return;
    }

    const QString name = QFileInfo(filePath).fileName();
    const QString base = name.toLower();
    const bool isPhimaxSnapshot = base.startsWith("phimax_") || base.startsWith("phimaxq4q6_") || base.startsWith("phimaxq");

    if (base == "run_config.txt") {
        QTextStream in(&f);
        struct Row {
            QString key;
            QString value;
        };
        QVector<Row> rows;
        rows.reserve(64);
        while (!in.atEnd()) {
            const QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;
            if (line.startsWith('#')) continue;
            const QStringList parts = line.split(QRegularExpression(R"(\s+)"), Qt::SkipEmptyParts);
            if (parts.isEmpty()) continue;
            const QString key = parts.value(0);
            const QString value = parts.mid(1).join(' ');
            rows.push_back(Row{key, value});
        }

        if (rows.isEmpty()) {
            resultsStatus_->setText("No config entries found in: " + name);
            return;
        }

        resultsTable_->setRowCount(rows.size());
        resultsTable_->setColumnCount(2);
        resultsTable_->setHorizontalHeaderLabels(QStringList() << "key" << "value");
        if (auto* header = resultsTable_->horizontalHeader()) {
            header->setStretchLastSection(false);
            header->setSectionResizeMode(QHeaderView::ResizeToContents);
        }
        for (int r = 0; r < rows.size(); ++r) {
            auto* keyItem = new QTableWidgetItem(rows[r].key);
            auto* valItem = new QTableWidgetItem(rows[r].value);
            valItem->setTextAlignment(Qt::AlignVCenter | Qt::AlignRight);
            resultsTable_->setItem(r, 0, keyItem);
            resultsTable_->setItem(r, 1, valItem);
        }
        resultsTable_->resizeColumnsToContents();
        resultsStatus_->setText(QString("%1 • %2 param(s)").arg(name).arg(rows.size()));
        return;
    }

    if (base.contains("checkpoint") && base.contains("timestamp")) {
        QTextStream in(&f);
        struct Entry {
            int step = 0;
            qint64 timestampMs = 0;
        };
        QVector<Entry> entries;
        entries.reserve(128);

        static const QRegularExpression re(R"(step\s+(\d+)\s+timestamp_ms\s+(\d+))",
                                           QRegularExpression::CaseInsensitiveOption);

        while (!in.atEnd()) {
            const QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;
            if (line.startsWith('#')) continue;

            auto m = re.match(line);
            if (m.hasMatch()) {
                bool ok1 = false;
                bool ok2 = false;
                const int step = m.captured(1).toInt(&ok1);
                const qint64 ts = m.captured(2).toLongLong(&ok2);
                if (ok1 && ok2) entries.push_back(Entry{step, ts});
                continue;
            }

            const QVector<double> nums = extractNumbersFromLine(line);
            if (nums.size() < 2) continue;
            const int step = static_cast<int>(nums[0]);
            const qint64 ts = static_cast<qint64>(nums[1]);
            entries.push_back(Entry{step, ts});
        }

        if (entries.isEmpty()) {
            resultsStatus_->setText("No checkpoint timestamp entries found in: " + name);
            return;
        }

        const qint64 t0 = entries.first().timestampMs;
        resultsTable_->setRowCount(entries.size());
        resultsTable_->setColumnCount(5);
        resultsTable_->setHorizontalHeaderLabels(
            QStringList() << "step"
                          << "timestamp_ms"
                          << "datetime"
                          << "elapsed_s"
                          << "delta_s");
        if (auto* header = resultsTable_->horizontalHeader()) {
            header->setStretchLastSection(false);
            header->setSectionResizeMode(QHeaderView::ResizeToContents);
        }

        qint64 prevTs = 0;
        for (int r = 0; r < entries.size(); ++r) {
            const auto& e = entries[r];
            const QDateTime dt = QDateTime::fromMSecsSinceEpoch(e.timestampMs);
            const QString dtStr = dt.toString("yyyy-MM-dd HH:mm:ss.zzz");
            const double elapsedS = static_cast<double>(e.timestampMs - t0) / 1000.0;
            const QString deltaStr =
                (r == 0 || prevTs <= 0) ? QString() : QString::number(static_cast<double>(e.timestampMs - prevTs) / 1000.0, 'f', 3);

            auto* stepItem = new QTableWidgetItem(QString::number(e.step));
            stepItem->setTextAlignment(Qt::AlignVCenter | Qt::AlignRight);
            auto* tsItem = new QTableWidgetItem(QString::number(e.timestampMs));
            tsItem->setTextAlignment(Qt::AlignVCenter | Qt::AlignRight);
            auto* dtItem = new QTableWidgetItem(dtStr);
            auto* elItem = new QTableWidgetItem(QString::number(elapsedS, 'f', 3));
            elItem->setTextAlignment(Qt::AlignVCenter | Qt::AlignRight);
            auto* dItem = new QTableWidgetItem(deltaStr);
            dItem->setTextAlignment(Qt::AlignVCenter | Qt::AlignRight);

            resultsTable_->setItem(r, 0, stepItem);
            resultsTable_->setItem(r, 1, tsItem);
            resultsTable_->setItem(r, 2, dtItem);
            resultsTable_->setItem(r, 3, elItem);
            resultsTable_->setItem(r, 4, dItem);

            prevTs = e.timestampMs;
        }
        resultsTable_->resizeColumnsToContents();
        const double totalS = static_cast<double>(entries.last().timestampMs - t0) / 1000.0;
        resultsStatus_->setText(QString("%1 • %2 checkpoint(s) • total %3 s").arg(name).arg(entries.size()).arg(totalS, 0, 'f', 3));
        return;
    }

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
        resultsStatus_->setText("No numeric data found in: " + name);
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

    QStringList headers;
    headers.reserve(maxCols);
    for (int c = 0; c < maxCols; ++c) headers << colLabel(c);
    resultsTable_->setHorizontalHeaderLabels(headers);

    if (isPhimaxSnapshot) {
        if (resultsImageValue_) {
            resultsImageValue_->clear();
            if (maxCols >= 4) resultsImageValue_->addItem("phimax", 3);
            if (maxCols >= 5) resultsImageValue_->addItem("con", 4);
            if (maxCols >= 6) resultsImageValue_->addItem("q4", 5);
            if (maxCols >= 7) resultsImageValue_->addItem("q6", 6);
            resultsImageValue_->setEnabled(resultsImageValue_->count() > 0);
        }
        if (resultsImagePointSize_) {
            resultsImagePointSize_->setEnabled(true);
        }
        resultsImageEnabled_ = resultsImageValue_ && (resultsImageValue_->count() > 0);
        if (resultsPlotPage_) {
            const int plotTab = resultsViewTabs_->indexOf(resultsPlotPage_);
            if (plotTab >= 0) {
                resultsViewTabs_->setTabEnabled(plotTab, true);
                if (resultsViewTabs_->tabBar()) resultsViewTabs_->tabBar()->setTabVisible(plotTab, true);
            }
        }
        resultsStatus_->setText(QString("%1 • %2 row(s) • %3 col(s) • use Plot tab").arg(name).arg(rows.size()).arg(maxCols));
        if (resultsPlotPage_) resultsViewTabs_->setCurrentWidget(resultsPlotPage_);
        renderResultsImage();
        return;
    }

    resultsStatus_->setText(QString("%1 • %2 row(s) • %3 col(s) • table only").arg(name).arg(rows.size()).arg(maxCols));
    if (resultsViewTabs_->count() > 0) {
        const int plotTab = resultsPlotPage_ ? resultsViewTabs_->indexOf(resultsPlotPage_) : -1;
        const int tableTab = (plotTab == 0) ? 1 : 0;
        if (tableTab >= 0 && tableTab < resultsViewTabs_->count()) resultsViewTabs_->setCurrentIndex(tableTab);
    }
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
        appendLog("=== Visualizer: Export requested but no image is available ===");
        QMessageBox::information(this, "Export image", "No image to export. Select a Phimax_*.txt file first.");
        return;
    }

    const QString baseDir = resultsDir_ ? resultsDir_->text().trimmed() : QString();
    const QString suggestDir = baseDir.isEmpty() ? QStandardPaths::writableLocation(QStandardPaths::PicturesLocation) : baseDir;
    const QString baseName = QFileInfo(resultsCurrentFile_).completeBaseName();
    const QString suggested = QDir(suggestDir).filePath(baseName + ".png");

    const QString fileName = pickSavePngFile(this, "Export PNG", suggested);
    if (fileName.isEmpty()) {
        appendLog("=== Visualizer: Export cancelled ===");
        return;
    }
    appendLog("=== Visualizer: Export saving to: " + fileName + " ===");

    if (!resultsImageItem_->pixmap().toImage().save(fileName)) {
        appendLog("=== Visualizer: Export failed ===");
        QMessageBox::warning(this, "Export image", "Failed to save: " + fileName);
        return;
    }
    appendLog("=== Visualizer: Export OK ===");
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

static QStringList splitValueList(const QString& text) {
    if (text.trimmed().isEmpty()) return {};
    QString s = text;
    s.replace(QChar(0xFF0C), QChar(',')); // Chinese comma
    s.replace(';', ',');
    const QStringList parts = s.split(QRegularExpression(R"([,\s]+)"), Qt::SkipEmptyParts);
    QStringList out;
    out.reserve(parts.size());
    for (const auto& p : parts) {
        const QString t = p.trimmed();
        if (!t.isEmpty()) out.push_back(t);
    }
    return out;
}

static QVector<double> parseDoubleList(const QString& text, QString* err) {
    QVector<double> out;
    const QStringList parts = splitValueList(text);
    if (parts.isEmpty()) {
        if (err) *err = "List is empty.";
        return out;
    }
    out.reserve(parts.size());
    for (const auto& s : parts) {
        bool ok = false;
        const double v = s.toDouble(&ok);
        if (!ok || !std::isfinite(v)) {
            if (err) *err = "Invalid number: " + s;
            out.clear();
            return out;
        }
        out.push_back(v);
    }
    return out;
}

static QVector<int> parseIntList(const QString& text, QString* err) {
    QVector<int> out;
    const QStringList parts = splitValueList(text);
    if (parts.isEmpty()) {
        if (err) *err = "List is empty.";
        return out;
    }
    out.reserve(parts.size());
    for (const auto& s : parts) {
        bool ok = false;
        const qlonglong v = s.toLongLong(&ok);
        if (!ok) {
            if (err) *err = "Invalid integer: " + s;
            out.clear();
            return out;
        }
        if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()) {
            if (err) *err = "Integer out of range: " + s;
            out.clear();
            return out;
        }
        out.push_back(static_cast<int>(v));
    }
    return out;
}

static QString joinPreviewDoubles(const QVector<double>& values, int maxItems = 64) {
    const int shown = std::min(static_cast<int>(values.size()), maxItems);
    QStringList parts;
    parts.reserve(shown);
    for (int i = 0; i < shown; ++i) parts << QString::number(values.at(i), 'g', 12);
    QString s = parts.join(',');
    if (values.size() > shown) s += QString(",…(+%1)").arg(values.size() - shown);
    return s;
}

static QString joinPreviewInts(const QVector<int>& values, int maxItems = 64) {
    const int shown = std::min(static_cast<int>(values.size()), maxItems);
    QStringList parts;
    parts.reserve(shown);
    for (int i = 0; i < shown; ++i) parts << QString::number(values.at(i));
    QString s = parts.join(',');
    if (values.size() > shown) s += QString(",…(+%1)").arg(values.size() - shown);
    return s;
}

QVector<RunJob> MainWindow::buildExperimentJobs(const QString& baseDir, const QString& experimentName, QString* errorOut) const {
    const int maxRuns = 5000;
    if (baseDir.trimmed().isEmpty()) {
        if (errorOut) *errorOut = "Base output dir is empty.";
        return {};
    }
    if (experimentName.trimmed().isEmpty()) {
        if (errorOut) *errorOut = "Experiment folder name is empty.";
        return {};
    }
    const QString experimentRoot = QDir(baseDir).filePath(experimentName);

    struct ParamValues {
        QString key;
        bool isInt = false;
        QVector<double> dvals;
        QVector<int> ivals;
    };

    QVector<ParamValues> params;
    params.reserve(sweepParams_.size());

    auto addValuesOrFail = [&](const SweepParamWidgets& p) -> bool {
        if (!p.mode || !p.stack) return true;

        const int mode = p.mode->currentIndex(); // 0 fixed, 1 range, 2 list
        QString err;
        ParamValues pv;
        pv.key = p.key;
        pv.isInt = p.isInt;

        if (p.isInt) {
            const int minV = p.fixedInt ? p.fixedInt->minimum() : std::numeric_limits<int>::min();
            const int maxV = p.fixedInt ? p.fixedInt->maximum() : std::numeric_limits<int>::max();
            if (mode == 0) {
                pv.ivals = QVector<int>{p.fixedInt ? p.fixedInt->value() : 0};
            } else if (mode == 1) {
                pv.ivals = makeIntRange(p.rangeStartInt->value(), p.rangeEndInt->value(), p.rangeStepInt->value(), &err);
            } else {
                pv.ivals = parseIntList(p.listEdit ? p.listEdit->text() : QString(), &err);
            }
            if (!err.isEmpty()) {
                if (errorOut) *errorOut = p.key + ": " + err;
                return false;
            }
            for (const int v : pv.ivals) {
                if (v < minV || v > maxV) {
                    if (errorOut) {
                        *errorOut = QString("%1: value %2 out of range [%3, %4].")
                                        .arg(p.key)
                                        .arg(v)
                                        .arg(minV)
                                        .arg(maxV);
                    }
                    return false;
                }
            }
        } else {
            const double minV = p.fixedDouble ? p.fixedDouble->minimum() : -std::numeric_limits<double>::infinity();
            const double maxV = p.fixedDouble ? p.fixedDouble->maximum() : std::numeric_limits<double>::infinity();
            if (mode == 0) {
                pv.dvals = QVector<double>{p.fixedDouble ? p.fixedDouble->value() : 0.0};
            } else if (mode == 1) {
                pv.dvals = makeDoubleRange(p.rangeStartDouble->value(), p.rangeEndDouble->value(), p.rangeStepDouble->value(), &err);
            } else {
                pv.dvals = parseDoubleList(p.listEdit ? p.listEdit->text() : QString(), &err);
            }
            if (!err.isEmpty()) {
                if (errorOut) *errorOut = p.key + ": " + err;
                return false;
            }
            for (const double v : pv.dvals) {
                if (!(v >= minV - 1e-12 && v <= maxV + 1e-12)) {
                    if (errorOut) {
                        *errorOut = QString("%1: value %2 out of range [%3, %4].")
                                        .arg(p.key)
                                        .arg(QString::number(v, 'g', 12))
                                        .arg(QString::number(minV, 'g', 12))
                                        .arg(QString::number(maxV, 'g', 12));
                    }
                    return false;
                }
            }
        }

        params.push_back(std::move(pv));
        return true;
    };

    for (const auto& p : sweepParams_) {
        if (!addValuesOrFail(p)) return {};
    }

    long long total = 1;
    for (const auto& p : params) {
        const int count = p.isInt ? p.ivals.size() : p.dvals.size();
        if (count <= 0) {
            if (errorOut) *errorOut = "No runs to execute.";
            return {};
        }
        if (total > (maxRuns / std::max(1, count))) {
            if (errorOut) *errorOut = QString("Too many runs (> %1). Reduce ranges/lists (limit: %1).").arg(maxRuns);
            return {};
        }
        total *= count;
    }
    if (total <= 0) {
        if (errorOut) *errorOut = "No runs to execute.";
        return {};
    }

    auto applyParamDouble = [](RunParams& rp, const QString& key, double v) {
        if (key == "u0") rp.u0 = v;
        else if (key == "con0") rp.con0 = v;
        else if (key == "sig") rp.sig = v;
        else if (key == "dt") rp.dt = v;
        else if (key == "dx") rp.dx = v;
    };
    auto applyParamInt = [](RunParams& rp, const QString& key, int v) {
        if (key == "steps") rp.steps = v;
        else if (key == "mod") rp.mod = v;
        else if (key == "seed") rp.seed = v;
    };

    QVector<int> indices;
    indices.fill(0, params.size());

    QVector<RunJob> jobs;
    jobs.reserve(static_cast<int>(total));

    for (long long i = 0; i < total; ++i) {
        RunParams rp;
        for (int pi = 0; pi < params.size(); ++pi) {
            const auto& pv = params.at(pi);
            if (pv.isInt) applyParamInt(rp, pv.key, pv.ivals.at(indices.at(pi)));
            else applyParamDouble(rp, pv.key, pv.dvals.at(indices.at(pi)));
        }

        const int idx = static_cast<int>(i) + 1;
        RunJob job;
        job.mode = "experiment";
        job.index = idx;
        job.total = static_cast<int>(total);
        job.params = rp;
        if (total <= 1) {
            job.outDir = experimentRoot;
        } else {
            job.outDir = QDir(experimentRoot).filePath(QString("run_%1").arg(idx, 4, 10, QChar('0')));
        }
        jobs.push_back(job);

        // Odometer increment (last param changes fastest).
        for (int k = params.size() - 1; k >= 0; --k) {
            const int count = params.at(k).isInt ? params.at(k).ivals.size() : params.at(k).dvals.size();
            indices[k]++;
            if (indices[k] < count) break;
            indices[k] = 0;
        }
    }

    return jobs;
}

void MainWindow::appendLog(const QString& text) {
    if (log_) log_->appendPlainText(text);

    QString single = text;
    single.replace('\n', ' ');
    single = single.simplified();
    statusLogText_ = single;
    statusLogTooltip_ = text;
    updateStatusLogText();
}

void MainWindow::updateStatusLogText() {
    if (!logStatusButton_) return;

    const QString full = statusLogText_.isEmpty() ? QString("Ready") : statusLogText_;
    int available = logStatusButton_->width();
    if (available <= 0 && statusBar()) available = statusBar()->width();
    available = std::max(60, available - 40);

    const QFontMetrics fm(logStatusButton_->font());
    logStatusButton_->setText(fm.elidedText(full, Qt::ElideRight, available));
    logStatusButton_->setToolTip(statusLogTooltip_.isEmpty() ? full : statusLogTooltip_);
}

void MainWindow::setRunningUi(bool running) {
    (void)running;
}

void MainWindow::resizeEvent(QResizeEvent* event) {
    QMainWindow::resizeEvent(event);
    updateStatusLogText();
}

void MainWindow::launchJob(const RunJob& job) {
    RunJob actual = job;

    QString err;
    if (!ensureDir(actual.outDir, &err)) {
        QMessageBox::critical(this, "Output error", err);
        return;
    }
    if (!writeParamsJson(actual, &err)) {
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
    process_->setArguments(buildArgs(actual.params, actual.outDir));
    process_->setWorkingDirectory(actual.outDir);
    process_->setProcessChannelMode(QProcess::MergedChannels);
    processOutputBuffer_.clear();
    currentJobMode_ = actual.mode;
    currentJobTotalSteps_ = std::max(0, actual.params.steps);
    if (expStepProgress_) {
        expStepProgress_->setRange(0, std::max(1, actual.params.steps));
        expStepProgress_->setValue(0);
    }

    connect(process_, &QProcess::readyRead, this, &MainWindow::onProcessReadyRead);
    connect(process_, &QProcess::finished, this, &MainWindow::onProcessFinished);
    connect(process_, &QProcess::errorOccurred, this, &MainWindow::onProcessError);

    currentOutputDir_ = actual.outDir;
    appendLog(QString("=== %1: %2 ===").arg(formatRunLabel(actual.index, actual.total), actual.outDir));
    appendLog(program + " " + process_->arguments().join(' '));

    process_->start();
}

void MainWindow::startExperiment() {
    if (process_) {
        QMessageBox::warning(this, "Busy", "A run is already in progress.");
        return;
    }
    const QString baseDir = expBaseOutDir_ ? expBaseOutDir_->text().trimmed() : QString();
    const QString experimentName =
        (expName_ && !expName_->text().trimmed().isEmpty()) ? expName_->text().trimmed() : makeTimestampedDirName("exp");

    QString err;
    jobQueue_ = buildExperimentJobs(baseDir, experimentName, &err);
    if (!err.isEmpty()) {
        QMessageBox::critical(this, "Experiment config error", err);
        return;
    }
    if (jobQueue_.isEmpty()) {
        QMessageBox::warning(this, "Experiment", "No runs to execute.");
        return;
    }
    currentJobIndex_ = 0;
    if (expProgress_) expProgress_->setValue(0);
    launchJob(jobQueue_.at(currentJobIndex_));
}

void MainWindow::stopRun() {
    jobQueue_.clear();
    currentJobIndex_ = -1;
    currentJobMode_.clear();
    currentJobTotalSteps_ = 0;
    if (expStepProgress_) expStepProgress_->setValue(0);
    if (expProgress_) expProgress_->setValue(0);
    if (!process_) return;
    appendLog("=== Stop requested ===");
    process_->kill();
}

void MainWindow::browseExperimentOutDir() {
    if (!expBaseOutDir_) return;
    const QString base = expBaseOutDir_->text().trimmed();
    const QString dir = pickExistingDirectory(this, "Select output directory", base.isEmpty() ? defaultBaseOutputDir() : base);
    if (!dir.isEmpty()) expBaseOutDir_->setText(dir);
}

void MainWindow::openCurrentOutputDir() {
    const QString dir = currentOutputDir_.isEmpty() ? (expBaseOutDir_ ? expBaseOutDir_->text() : QString()) : currentOutputDir_;
    if (dir.trimmed().isEmpty()) return;
    QDesktopServices::openUrl(QUrl::fromLocalFile(dir));
}

void MainWindow::updateExperimentPreview() {
    if (!expPreview_) return;
    const int maxRuns = 5000;

    long long total = 1;
    QStringList varied;

    for (auto& p : sweepParams_) {
        if (!p.mode || !p.stack) continue;

        const int mode = p.mode->currentIndex(); // 0 fixed, 1 range, 2 list
        p.stack->setCurrentIndex(mode);

        QString err;
        int count = 1;

        if (p.isInt) {
            QVector<int> vals;
            const int minV = p.fixedInt ? p.fixedInt->minimum() : std::numeric_limits<int>::min();
            const int maxV = p.fixedInt ? p.fixedInt->maximum() : std::numeric_limits<int>::max();

            if (mode == 0) {
                vals = QVector<int>{p.fixedInt ? p.fixedInt->value() : 0};
            } else if (mode == 1) {
                vals = makeIntRange(p.rangeStartInt->value(), p.rangeEndInt->value(), p.rangeStepInt->value(), &err);
                if (p.rangePreview) p.rangePreview->setText(joinPreviewInts(vals));
            } else {
                if (p.rangePreview) p.rangePreview->clear();
                vals = parseIntList(p.listEdit ? p.listEdit->text() : QString(), &err);
            }

            if (!err.isEmpty()) {
                expPreview_->setText("Preview: " + p.key + ": " + err);
                return;
            }
            for (const int v : vals) {
                if (v < minV || v > maxV) {
                    expPreview_->setText(QString("Preview: %1: value %2 out of range [%3, %4].")
                                               .arg(p.key)
                                               .arg(v)
                                               .arg(minV)
                                               .arg(maxV));
                    return;
                }
            }
            count = vals.size();
        } else {
            QVector<double> vals;
            const double minV = p.fixedDouble ? p.fixedDouble->minimum() : -std::numeric_limits<double>::infinity();
            const double maxV = p.fixedDouble ? p.fixedDouble->maximum() : std::numeric_limits<double>::infinity();

            if (mode == 0) {
                vals = QVector<double>{p.fixedDouble ? p.fixedDouble->value() : 0.0};
            } else if (mode == 1) {
                vals = makeDoubleRange(p.rangeStartDouble->value(), p.rangeEndDouble->value(), p.rangeStepDouble->value(), &err);
                if (p.rangePreview) p.rangePreview->setText(joinPreviewDoubles(vals));
            } else {
                if (p.rangePreview) p.rangePreview->clear();
                vals = parseDoubleList(p.listEdit ? p.listEdit->text() : QString(), &err);
            }

            if (!err.isEmpty()) {
                expPreview_->setText("Preview: " + p.key + ": " + err);
                return;
            }
            for (const double v : vals) {
                if (!(v >= minV - 1e-12 && v <= maxV + 1e-12)) {
                    expPreview_->setText(QString("Preview: %1: value %2 out of range [%3, %4].")
                                               .arg(p.key)
                                               .arg(QString::number(v, 'g', 12))
                                               .arg(QString::number(minV, 'g', 12))
                                               .arg(QString::number(maxV, 'g', 12)));
                    return;
                }
            }
            count = vals.size();
        }

        if (count <= 0) {
            expPreview_->setText("Preview: No runs to execute.");
            return;
        }

        if (mode != 0 && count > 1) varied << QString("%1=%2").arg(p.key).arg(count);

        if (total <= maxRuns && total > (maxRuns / std::max(1, count))) {
            total = maxRuns + 1;
        } else if (total <= maxRuns) {
            total *= count;
        }
    }

    QString msg;
    if (total > maxRuns) {
        msg = QString("Preview: too many runs (> %1). Reduce ranges/lists.").arg(maxRuns);
    } else {
        msg = QString("Preview: %1 run(s).").arg(total);
    }
    if (!varied.isEmpty()) msg += " " + varied.join(", ");
    expPreview_->setText(msg);
}

void MainWindow::showAboutDialog() {
    const QString appName = QCoreApplication::applicationName().isEmpty() ? QString("MID Nano") : QCoreApplication::applicationName();
    const QString version = QCoreApplication::applicationVersion().isEmpty() ? QString("dev") : QCoreApplication::applicationVersion();

    static const QString githubUrl = "https://github.com/zhiqing0205/PFC-Exp-GUI";
    static const QString giteeUrl = "https://gitee.com/zhiqing0205/PFC-Exp-GUI";
    static const QString author = "zhiqing0205";

    QDialog dlg(this);
    dlg.setWindowTitle("About " + appName);
    dlg.setModal(true);
    dlg.resize(560, 520);

    auto* layout = new QVBoxLayout(&dlg);
    layout->setSpacing(12);

    auto* title = new QLabel(QString("<div style='font-size:18pt; font-weight:700;'>%1</div>"
                                     "<div style='font-size:11pt; color:#6B7280;'>Version %2</div>")
                                 .arg(appName.toHtmlEscaped())
                                 .arg(version.toHtmlEscaped()));
    title->setTextFormat(Qt::RichText);
    title->setAlignment(Qt::AlignHCenter);
    layout->addWidget(title);

    auto makeKeyLabel = [](const QString& text) {
        auto* l = new QLabel(text);
        QFont f = l->font();
        f.setBold(true);
        l->setFont(f);
        l->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        return l;
    };
    auto makeValueLabel = [](const QString& html) {
        auto* l = new QLabel(html);
        l->setTextFormat(Qt::RichText);
        l->setOpenExternalLinks(true);
        l->setWordWrap(true);
        l->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::LinksAccessibleByMouse);
        return l;
    };

    auto makeGroup = [&](const QString& titleText) {
        auto* box = new QGroupBox(titleText);
        auto* form = new QFormLayout;
        form->setLabelAlignment(Qt::AlignRight);
        form->setFormAlignment(Qt::AlignTop);
        box->setLayout(form);
        return std::pair<QGroupBox*, QFormLayout*>(box, form);
    };

    auto* scroll = new QScrollArea;
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setWidgetResizable(true);

    auto* content = new QWidget;
    auto* contentLayout = new QVBoxLayout(content);
    contentLayout->setContentsMargins(0, 0, 0, 0);
    contentLayout->setSpacing(12);

    // License
    {
        const QString appDataDir = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
        bool installed = false;
        if (!appDataDir.isEmpty()) {
            QDir dir(appDataDir);
            const QStringList files = dir.entryList(QStringList() << "license.*", QDir::Files);
            installed = !files.isEmpty();
        }
        const QString statusHtml = installed ? QString("<span style='color:#16A34A; font-weight:700;'>Installed</span>")
                                             : QString("<span style='color:#DC2626; font-weight:700;'>Not installed</span>");
        auto [box, f] = makeGroup("License");
        f->addRow(makeKeyLabel("Status"), makeValueLabel(statusHtml));
        contentLayout->addWidget(box);
    }

    // Project
    {
        auto [box, f] = makeGroup("Project");
        f->addRow(makeKeyLabel("Author"), makeValueLabel(author.toHtmlEscaped()));
        contentLayout->addWidget(box);
    }

    // Links
    {
        auto [box, f] = makeGroup("Links");
        f->addRow(makeKeyLabel("GitHub"), makeValueLabel(QString("<a href=\"%1\">%1</a>").arg(githubUrl.toHtmlEscaped())));
        f->addRow(makeKeyLabel("Gitee"), makeValueLabel(QString("<a href=\"%1\">%1</a>").arg(giteeUrl.toHtmlEscaped())));
        contentLayout->addWidget(box);
    }

    // Dependencies
    {
        auto [box, f] = makeGroup("Dependencies");
        f->addRow(makeKeyLabel("Qt"), makeValueLabel("<a href=\"https://www.qt.io\">https://www.qt.io</a>"));
        f->addRow(makeKeyLabel("FFTW"), makeValueLabel("<a href=\"https://www.fftw.org\">https://www.fftw.org</a>"));
        contentLayout->addWidget(box);
    }

    contentLayout->addStretch(1);
    scroll->setWidget(content);
    layout->addWidget(scroll, 1);

    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok);
    connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    layout->addWidget(buttons);

    dlg.exec();
}

void MainWindow::uploadLicense() {
    const QString selected = QFileDialog::getOpenFileName(
        this,
        "Select license file",
        QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation),
        "License files (*.lic *.txt *.json);;All files (*.*)");

    if (selected.trimmed().isEmpty()) {
        appendLog("=== License: upload cancelled ===");
        return;
    }

    const QString appDataDir = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
    if (appDataDir.isEmpty()) {
        QMessageBox::warning(this, "License", "Cannot determine an app data directory to store the license file.");
        return;
    }

    QDir dir(appDataDir);
    if (!dir.exists() && !dir.mkpath(".")) {
        QMessageBox::warning(this, "License", "Failed to create app data directory:\n" + appDataDir);
        return;
    }

    const QString suffix = QFileInfo(selected).suffix();
    const QString destName = suffix.isEmpty() ? QString("license.lic") : QString("license.%1").arg(suffix);
    const QString destPath = dir.filePath(destName);
    if (QFileInfo::exists(destPath)) QFile::remove(destPath);

    if (!QFile::copy(selected, destPath)) {
        QMessageBox::warning(this, "License", "Failed to copy license file to:\n" + destPath);
        appendLog("=== License: upload failed: " + destPath + " ===");
        return;
    }

    appendLog("=== License: saved to: " + destPath + " ===");
    QMessageBox::information(this, "License", "License file saved to:\n" + destPath);
}

void MainWindow::onProcessReadyRead() {
    if (!process_) return;
    const QByteArray data = process_->readAll();
    if (data.isEmpty()) return;

    processOutputBuffer_ += QString::fromLocal8Bit(data);

    static const QRegularExpression progressRe(R"(^PFC_PROGRESS\s+step=(\d+)\s+total=(\d+)\s*$)");
    int nl = -1;
    while ((nl = processOutputBuffer_.indexOf('\n')) >= 0) {
        QString line = processOutputBuffer_.left(nl);
        processOutputBuffer_.remove(0, nl + 1);
        if (line.endsWith('\r')) line.chop(1);

        const auto m = progressRe.match(line);
        if (m.hasMatch()) {
            const int step = m.captured(1).toInt();
            const int total = m.captured(2).toInt();

            if (expStepProgress_) {
                if (total > 0 && expStepProgress_->maximum() != total) {
                    expStepProgress_->setRange(0, total);
                }
                expStepProgress_->setValue(std::clamp(step, 0, expStepProgress_->maximum()));
            }
            continue;
        }

        if (!line.isEmpty()) appendLog(line);
    }
}

void MainWindow::onProcessFinished(int exitCode, QProcess::ExitStatus exitStatus) {
    onProcessReadyRead();
    if (!processOutputBuffer_.trimmed().isEmpty()) {
        appendLog(processOutputBuffer_);
    }
    processOutputBuffer_.clear();

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
        appendLog("=== Hint: SIGKILL(9) on macOS can be Gatekeeper/quarantine or a missing dylib. If this is a downloaded app, try: xattr -dr com.apple.quarantine \"/Applications/MID Nano.app\" ===");
    }
#endif

    process_->deleteLater();
    process_ = nullptr;

    if (ok) {
        if (expStepProgress_) expStepProgress_->setValue(expStepProgress_->maximum());
    }

    if (jobQueue_.isEmpty() || currentJobIndex_ < 0) {
        if (expProgress_) expProgress_->setValue(0);
        currentJobIndex_ = -1;
        currentJobMode_.clear();
        currentJobTotalSteps_ = 0;
        return;
    }

    const int total = jobQueue_.size();
    const int done = currentJobIndex_ + 1;
    if (expProgress_) expProgress_->setValue(static_cast<int>(100.0 * done / total));

    if (!ok) {
        appendLog("=== Experiment aborted due to failure ===");
        jobQueue_.clear();
        currentJobIndex_ = -1;
        return;
    }

    ++currentJobIndex_;
    if (currentJobIndex_ >= total) {
        appendLog("=== Experiment complete ===");
        jobQueue_.clear();
        currentJobIndex_ = -1;
        if (expProgress_) expProgress_->setValue(100);
        currentJobMode_.clear();
        currentJobTotalSteps_ = 0;
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

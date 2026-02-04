#include "MainWindow.h"

#include <QApplication>
#include <QIcon>
#include <QLinearGradient>
#include <QPainter>
#include <QPalette>
#include <QPixmap>
#include <QStyleFactory>

static QIcon makeAppIcon() {
    QIcon icon;
    const QList<int> sizes{16, 24, 32, 48, 64, 128, 256};
    for (const int size : sizes) {
        QPixmap pm(size, size);
        pm.fill(Qt::transparent);

        QPainter p(&pm);
        p.setRenderHint(QPainter::Antialiasing, true);

        const QRectF r(0, 0, size, size);
        QLinearGradient g(r.topLeft(), r.bottomRight());
        g.setColorAt(0.0, QColor(64, 128, 255));
        g.setColorAt(1.0, QColor(128, 80, 255));

        p.setPen(Qt::NoPen);
        p.setBrush(g);
        p.drawRoundedRect(r.adjusted(size * 0.08, size * 0.08, -size * 0.08, -size * 0.08), size * 0.22, size * 0.22);

        QFont f = p.font();
        f.setBold(true);
        f.setPixelSize(static_cast<int>(size * 0.30));
        p.setFont(f);
        p.setPen(Qt::white);
        p.drawText(r, Qt::AlignCenter, "MID\nNano");

        p.end();
        icon.addPixmap(pm);
    }
    return icon;
}

static QString appStyleSheet() {
    return QString::fromUtf8(R"QSS(
QWidget {
  color: #111827;
  font-size: 12pt;
}

QMainWindow {
  background: #F5F7FB;
}

QTabWidget::pane {
  border: 0;
  margin-top: -1px;
}

QTabBar::base {
  border: 0;
  background: #F5F7FB;
}

QTabWidget#mainTabs::right-corner, QTabWidget#mainTabs::left-corner {
  background: #F5F7FB;
  border: 0px;
}
QWidget#mainTabsActions {
  background: #F5F7FB;
  border: 1px solid #D8DEE9;
  border-radius: 12px;
}

QTabBar::tab {
  background: #EEF2F7;
  border: 1px solid #D8DEE9;
  border-bottom: 0;
  padding: 8px 14px;
  margin-right: 6px;
  border-top-left-radius: 10px;
  border-top-right-radius: 10px;
}
QTabBar::tab:selected {
  background: #FFFFFF;
  border-color: #CBD5E1;
}
QTabBar::tab:hover {
  background: #F4F6FB;
}

QGroupBox {
  border: 1px solid #D8DEE9;
  border-radius: 14px;
  margin-top: 18px;
  padding: 12px;
  background: #FFFFFF;
}
QGroupBox::title {
  subcontrol-origin: margin;
  left: 14px;
  padding: 0 8px;
  color: #1D4ED8;
  font-weight: 600;
}

QLabel[hint="true"] {
  background: #FFFFFF;
  border: 1px solid #D8DEE9;
  border-radius: 12px;
  padding: 10px 12px;
  color: #374151;
}

QScrollArea {
  background: transparent;
}

QSplitter::handle {
  background: #E5E7EB;
}
QSplitter::handle:hover {
  background: #D1D5DB;
}

QScrollBar:vertical {
  background: transparent;
  width: 10px;
  margin: 0px;
}
QScrollBar::handle:vertical {
  background: #CBD5E1;
  border-radius: 5px;
  min-height: 30px;
}
QScrollBar::handle:vertical:hover {
  background: #94A3B8;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
  height: 0px;
}
QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
  background: transparent;
}

QScrollBar:horizontal {
  background: transparent;
  height: 10px;
  margin: 0px;
}
QScrollBar::handle:horizontal {
  background: #CBD5E1;
  border-radius: 5px;
  min-width: 30px;
}
QScrollBar::handle:horizontal:hover {
  background: #94A3B8;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
  width: 0px;
}
QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {
  background: transparent;
}

QCheckBox {
  spacing: 8px;
}
QCheckBox::indicator {
  width: 18px;
  height: 18px;
  border-radius: 5px;
  background: #FFFFFF;
  border: 1px solid #CBD5E1;
}
QCheckBox::indicator:hover {
  border: 1px solid #2563EB;
}
QCheckBox::indicator:checked {
  background: #2563EB;
  border: 1px solid #2563EB;
  image: url(:/icons/check.svg);
}
QCheckBox::indicator:checked:hover {
  background: #3B82F6;
  border: 1px solid #3B82F6;
  image: url(:/icons/check.svg);
}

QLineEdit, QSpinBox, QDoubleSpinBox {
  background: #FFFFFF;
  border: 1px solid #CBD5E1;
  border-radius: 10px;
  selection-background-color: #2563EB;
}
QLineEdit {
  padding: 6px 8px;
}
QSpinBox, QDoubleSpinBox {
  padding: 6px 34px 6px 10px;
}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus {
  border: 1px solid #2563EB;
}

QSpinBox::up-button, QDoubleSpinBox::up-button {
  subcontrol-origin: border;
  subcontrol-position: top right;
  width: 26px;
  border-left: 1px solid #E5E7EB;
  border-top-right-radius: 10px;
  background: transparent;
}
QSpinBox::down-button, QDoubleSpinBox::down-button {
  subcontrol-origin: border;
  subcontrol-position: bottom right;
  width: 26px;
  border-left: 1px solid #E5E7EB;
  border-bottom-right-radius: 10px;
  background: transparent;
}
QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover,
QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover {
  background: #E5E7EB;
}
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {
  image: url(:/icons/chevron-up.svg);
  width: 12px;
  height: 12px;
}
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {
  image: url(:/icons/chevron-down.svg);
  width: 12px;
  height: 12px;
}

QPlainTextEdit {
  background: #FFFFFF;
  border: 1px solid #D8DEE9;
  border-radius: 14px;
  padding: 10px;
  selection-background-color: #2563EB;
}

QComboBox {
  background: #FFFFFF;
  border: 1px solid #CBD5E1;
  border-radius: 10px;
  padding: 6px 34px 6px 10px;
}
QComboBox:focus {
  border: 1px solid #2563EB;
}
QComboBox:on {
  background: #FFFFFF;
  border: 1px solid #2563EB;
}
QComboBox::drop-down {
  subcontrol-origin: padding;
  subcontrol-position: top right;
  width: 30px;
  border-left: 1px solid #E5E7EB;
}
QComboBox::down-arrow {
  image: url(:/icons/chevron-down.svg);
  width: 12px;
  height: 12px;
}
QComboBox[kind="mode"] {
  background: #F0FDFA;
  border: 1px solid #99F6E4;
}
QComboBox[kind="mode"]::drop-down {
  border-left: 1px solid #99F6E4;
}
QComboBox[kind="mode"]:hover {
  border-color: #14B8A6;
}
QComboBox[kind="mode"]:focus {
  border-color: #14B8A6;
}
QComboBox[kind="mode"]:on {
  background: #F0FDFA;
  border-color: #14B8A6;
}
QComboBox[kind="mode"]::drop-down:on {
  border-left: 1px solid #14B8A6;
}
QComboBox QAbstractItemView {
  background: #FFFFFF;
  border: 1px solid #CBD5E1;
  selection-background-color: #DBEAFE;
  selection-color: #111827;
}
QComboBox QAbstractItemView::item {
  padding: 6px 10px;
}
QComboBox[kind="mode"] QAbstractItemView {
  background: #F0FDFA;
  border: 1px solid #99F6E4;
  selection-background-color: #14B8A6;
  selection-color: #0B0B0B;
}

QToolButton, QPushButton {
  background: #EEF2F7;
  border: 1px solid #D1D5DB;
  border-radius: 10px;
  padding: 8px 14px;
}
QToolButton:hover, QPushButton:hover {
  background: #E5E7EB;
  border-color: #CBD5E1;
}
QToolButton:pressed, QPushButton:pressed {
  background: #D1D5DB;
}

QWidget#mainTabsActions QPushButton#tabUploadLicense,
QWidget#mainTabsActions QPushButton#tabAbout {
  background: #FFFFFF;
  border: 1px solid #CBD5E1;
  border-radius: 10px;
  padding: 6px 12px;
}
QWidget#mainTabsActions QPushButton#tabUploadLicense:hover,
QWidget#mainTabsActions QPushButton#tabAbout:hover {
  background: #E5E7EB;
  border-color: #CBD5E1;
}
QWidget#mainTabsActions QPushButton#tabUploadLicense:pressed,
QWidget#mainTabsActions QPushButton#tabAbout:pressed {
  background: #D1D5DB;
  border-color: #CBD5E1;
}

QPushButton[primary="true"] {
  background: #2563EB;
  border-color: #2563EB;
  color: #FFFFFF;
  font-weight: 700;
}
QPushButton[primary="true"]:hover {
  background: #3B82F6;
  border-color: #3B82F6;
}
QPushButton[danger="true"] {
  background: #DC2626;
  border-color: #DC2626;
  color: #FFFFFF;
  font-weight: 700;
}
QPushButton[danger="true"]:hover {
  background: #EF4444;
  border-color: #EF4444;
}

QProgressBar {
  background: #FFFFFF;
  border: 1px solid #D1D5DB;
  border-radius: 10px;
  padding: 2px;
  text-align: center;
}
QProgressBar::chunk {
  border-radius: 8px;
  background: #2563EB;
}

QListWidget {
  background: #FFFFFF;
  border: 1px solid #D8DEE9;
  border-radius: 14px;
  padding: 6px;
}
QListWidget::item {
  padding: 6px 8px;
  border-radius: 10px;
}
QListWidget::item:selected {
  background: #2563EB;
  color: #FFFFFF;
}
QListWidget::item:selected:hover {
  background: #1D4ED8;
}

QTableView {
  background: #FFFFFF;
  border: 1px solid #D8DEE9;
  border-radius: 14px;
  gridline-color: #E5E7EB;
  selection-background-color: #DBEAFE;
}
QHeaderView::section {
  background: #EEF2F7;
  color: #1D4ED8;
  padding: 6px 8px;
  border: 0px;
}

QStatusBar {
  background: #EEF2F7;
  border-top: 1px solid #D8DEE9;
}
QStatusBar::item {
  border: 0px;
}
QStatusBar QToolButton {
  background: transparent;
  border: 0px;
  padding: 6px 10px;
  border-radius: 10px;
}
QStatusBar QToolButton:hover {
  background: #E5E7EB;
}
)QSS");
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    app.setApplicationName("MID Nano");
#ifdef MID_NANO_VERSION
    app.setApplicationVersion(QString::fromUtf8(MID_NANO_VERSION));
#endif
    app.setWindowIcon(makeAppIcon());

    if (QStyleFactory::keys().contains("Fusion")) {
        QApplication::setStyle(QStyleFactory::create("Fusion"));
    }

    // A simple, modern light palette.
    QPalette palette;
    palette.setColor(QPalette::Window, QColor(245, 247, 251));
    palette.setColor(QPalette::WindowText, QColor(17, 24, 39));
    palette.setColor(QPalette::Base, QColor(255, 255, 255));
    palette.setColor(QPalette::AlternateBase, QColor(238, 242, 247));
    palette.setColor(QPalette::ToolTipBase, QColor(255, 255, 255));
    palette.setColor(QPalette::ToolTipText, QColor(17, 24, 39));
    palette.setColor(QPalette::Text, QColor(17, 24, 39));
    palette.setColor(QPalette::Button, QColor(238, 242, 247));
    palette.setColor(QPalette::ButtonText, QColor(17, 24, 39));
    palette.setColor(QPalette::BrightText, Qt::red);
    palette.setColor(QPalette::Highlight, QColor(37, 99, 235));
    palette.setColor(QPalette::HighlightedText, Qt::white);
    app.setPalette(palette);

    app.setStyleSheet(appStyleSheet());

    MainWindow w;
    w.show();
    return app.exec();
}

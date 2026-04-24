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

/* ── Tab bar: modern underline style ── */
QTabWidget::pane {
  border: 0;
  border-top: 1px solid #E5E7EB;
}
QTabWidget#mainTabs::pane {
  background: #F5F7FB;
}
QTabBar::tab {
  background: transparent;
  border: none;
  border-bottom: 3px solid transparent;
  padding: 10px 20px;
  margin-right: 2px;
  color: #6B7280;
  font-weight: 500;
}
QTabBar::tab:selected {
  color: #2563EB;
  border-bottom: 3px solid #2563EB;
  font-weight: 600;
}
QTabBar::tab:hover:!selected {
  color: #374151;
  border-bottom: 3px solid #D1D5DB;
}

/* ── Cards & GroupBoxes ── */
QGroupBox {
  border: 1px solid #E5E7EB;
  border-radius: 12px;
  margin-top: 20px;
  padding: 16px;
  background: #FFFFFF;
}
QGroupBox::title {
  subcontrol-origin: margin;
  left: 16px;
  padding: 0 10px;
  color: #1D4ED8;
  font-weight: 600;
  font-size: 13pt;
}

/* ── Hint labels: accent left border ── */
QLabel[hint="true"] {
  background: #F8FAFC;
  border: none;
  border-left: 3px solid #2563EB;
  border-radius: 0px;
  padding: 10px 14px;
  color: #374151;
  font-size: 11pt;
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

/* ── Scrollbars ── */
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

/* ── Checkboxes ── */
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

/* ── Input fields ── */
QLineEdit, QSpinBox, QDoubleSpinBox {
  background: #FFFFFF;
  border: 1px solid #CBD5E1;
  border-radius: 8px;
  selection-background-color: #2563EB;
}
QLineEdit {
  padding: 6px 10px;
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
  border-top-right-radius: 8px;
  background: transparent;
}
QSpinBox::down-button, QDoubleSpinBox::down-button {
  subcontrol-origin: border;
  subcontrol-position: bottom right;
  width: 26px;
  border-left: 1px solid #E5E7EB;
  border-bottom-right-radius: 8px;
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
  border: 1px solid #E5E7EB;
  border-radius: 12px;
  padding: 10px;
  selection-background-color: #2563EB;
}

/* ── ComboBox ── */
QComboBox {
  background: #FFFFFF;
  border: 1px solid #CBD5E1;
  border-radius: 8px;
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

/* ── Buttons ── */
QToolButton, QPushButton {
  background: #EEF2F7;
  border: 1px solid #D1D5DB;
  border-radius: 8px;
  padding: 10px 20px;
  font-weight: 500;
}
QToolButton:hover, QPushButton:hover {
  background: #E2E6ED;
  border-color: #B0B8C4;
}
QToolButton:pressed, QPushButton:pressed {
  background: #D1D5DB;
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

/* ── Progress bars ── */
QProgressBar {
  background: #F1F5F9;
  border: 1px solid #E2E8F0;
  border-radius: 8px;
  padding: 1px;
  text-align: center;
  font-size: 10pt;
  color: #475569;
  min-height: 22px;
}
QProgressBar::chunk {
  border-radius: 7px;
  background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
    stop:0 #2563EB, stop:1 #3B82F6);
}

/* ── Lists & Tables ── */
QListWidget {
  background: #FFFFFF;
  border: 1px solid #E5E7EB;
  border-radius: 12px;
  padding: 6px;
}
QListWidget::item {
  padding: 6px 8px;
  border-radius: 8px;
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
  border: 1px solid #E5E7EB;
  border-radius: 12px;
  gridline-color: #E5E7EB;
  selection-background-color: #DBEAFE;
}
QHeaderView::section {
  background: #F1F5F9;
  color: #1D4ED8;
  padding: 6px 10px;
  border: 0px;
  font-weight: 500;
}

/* ── Status bar ── */
QStatusBar {
  background: #F1F5F9;
  border-top: 1px solid #E5E7EB;
}
QStatusBar::item {
  border: 0px;
}
QStatusBar QToolButton {
  background: transparent;
  border: 0px;
  padding: 6px 10px;
  border-radius: 8px;
}
QStatusBar QToolButton:hover {
  background: #E5E7EB;
}

/* ── Model toggle segmented control ── */
QPushButton[modelToggle="true"] {
  border: 1px solid #D1D5DB;
  padding: 8px 24px;
  background: #F8F9FA;
  font-weight: 500;
  font-size: 13px;
}
QPushButton[modelToggle="true"][togglePos="left"] {
  border-top-left-radius: 8px;
  border-bottom-left-radius: 8px;
  border-top-right-radius: 0px;
  border-bottom-right-radius: 0px;
  margin-right: -1px;
}
QPushButton[modelToggle="true"][togglePos="right"] {
  border-top-left-radius: 0px;
  border-bottom-left-radius: 0px;
  border-top-right-radius: 8px;
  border-bottom-right-radius: 8px;
}
QPushButton[modelToggle="true"]:checked {
  background: #2563EB;
  color: white;
  border-color: #2563EB;
}
QPushButton[modelToggle="true"]:hover:!checked {
  background: #E5E7EB;
  border-color: #9CA3AF;
}
)QSS");
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    app.setApplicationName("MID Nano");
    app.setOrganizationName("MIDNano");
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

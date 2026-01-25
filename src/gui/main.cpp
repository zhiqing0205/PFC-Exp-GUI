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
        f.setPixelSize(static_cast<int>(size * 0.42));
        p.setFont(f);
        p.setPen(Qt::white);
        p.drawText(r, Qt::AlignCenter, "PFC");

        p.end();
        icon.addPixmap(pm);
    }
    return icon;
}

static QString appStyleSheet() {
    return QString::fromUtf8(R"QSS(
QWidget {
  color: #EAEAEA;
  font-size: 11pt;
}

QMainWindow {
  background: #121212;
}

QTabWidget::pane {
  border: 0;
  margin-top: -1px;
}

QTabBar::base {
  border: 0;
  background: transparent;
}

QTabBar::tab {
  background: #1A1A1A;
  border: 1px solid #2A2A2A;
  border-bottom: 0;
  padding: 8px 14px;
  margin-right: 6px;
  border-top-left-radius: 10px;
  border-top-right-radius: 10px;
}
QTabBar::tab:selected {
  background: #202020;
  border-color: #3A3A3A;
}
QTabBar::tab:hover {
  background: #242424;
}

QGroupBox {
  border: 1px solid #2A2A2A;
  border-radius: 14px;
  margin-top: 18px;
  padding: 12px;
  background: #161616;
}
QGroupBox::title {
  subcontrol-origin: margin;
  left: 14px;
  padding: 0 8px;
  color: #B7C6FF;
  font-weight: 600;
}

QLabel[hint="true"] {
  background: #161616;
  border: 1px solid #2A2A2A;
  border-radius: 12px;
  padding: 10px 12px;
  color: #CFCFCF;
}

QScrollArea {
  background: transparent;
}

QSplitter::handle {
  background: #1A1A1A;
}
QSplitter::handle:hover {
  background: #242424;
}

QScrollBar:vertical {
  background: transparent;
  width: 10px;
  margin: 0px;
}
QScrollBar::handle:vertical {
  background: #2A2A2A;
  border-radius: 5px;
  min-height: 30px;
}
QScrollBar::handle:vertical:hover {
  background: #3A3A3A;
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
  background: #2A2A2A;
  border-radius: 5px;
  min-width: 30px;
}
QScrollBar::handle:horizontal:hover {
  background: #3A3A3A;
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
  background: #101010;
  border: 1px solid #2E2E2E;
}
QCheckBox::indicator:hover {
  border: 1px solid #4080FF;
}
QCheckBox::indicator:checked {
  background: #4080FF;
  border: 1px solid #4080FF;
}
QCheckBox::indicator:checked:hover {
  background: #5B93FF;
  border: 1px solid #5B93FF;
}

QLineEdit, QSpinBox, QDoubleSpinBox {
  background: #101010;
  border: 1px solid #2E2E2E;
  border-radius: 10px;
  padding: 6px 8px;
  selection-background-color: #4080FF;
}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus {
  border: 1px solid #4080FF;
}

QPlainTextEdit {
  background: #0E0E0E;
  border: 1px solid #2A2A2A;
  border-radius: 14px;
  padding: 10px;
  selection-background-color: #4080FF;
}

QComboBox {
  background: #101010;
  border: 1px solid #2E2E2E;
  border-radius: 10px;
  padding: 6px 10px;
}
QComboBox:focus {
  border: 1px solid #4080FF;
}
QComboBox::drop-down {
  border: 0px;
  width: 26px;
}
QComboBox QAbstractItemView {
  background: #101010;
  border: 1px solid #2E2E2E;
  selection-background-color: #4080FF;
}

QToolButton, QPushButton {
  background: #232323;
  border: 1px solid #323232;
  border-radius: 10px;
  padding: 8px 14px;
}
QToolButton:hover, QPushButton:hover {
  background: #2A2A2A;
  border-color: #3A3A3A;
}
QToolButton:pressed, QPushButton:pressed {
  background: #1B1B1B;
}

QPushButton[primary="true"] {
  background: #4080FF;
  border-color: #4080FF;
  color: #0B0B0B;
  font-weight: 700;
}
QPushButton[primary="true"]:hover {
  background: #5B93FF;
  border-color: #5B93FF;
}
QPushButton[danger="true"] {
  background: #E05252;
  border-color: #E05252;
  color: #0B0B0B;
  font-weight: 700;
}
QPushButton[danger="true"]:hover {
  background: #F06A6A;
  border-color: #F06A6A;
}

QProgressBar {
  background: #101010;
  border: 1px solid #2E2E2E;
  border-radius: 10px;
  padding: 2px;
  text-align: center;
}
QProgressBar::chunk {
  border-radius: 8px;
  background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
    stop:0 #4080FF, stop:1 #7A5CFF);
}

QListWidget {
  background: #0E0E0E;
  border: 1px solid #2A2A2A;
  border-radius: 14px;
  padding: 6px;
}
QListWidget::item {
  padding: 6px 8px;
  border-radius: 10px;
}
QListWidget::item:selected {
  background: #253A68;
}

QTableView {
  background: #0E0E0E;
  border: 1px solid #2A2A2A;
  border-radius: 14px;
  gridline-color: #2A2A2A;
  selection-background-color: #253A68;
}
QHeaderView::section {
  background: #1A1A1A;
  color: #B7C6FF;
  padding: 6px 8px;
  border: 0px;
}
)QSS");
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    app.setApplicationName("PFC-Exp-GUI");
    app.setWindowIcon(makeAppIcon());

    if (QStyleFactory::keys().contains("Fusion")) {
        QApplication::setStyle(QStyleFactory::create("Fusion"));
    }

    // A simple, modern dark palette.
    QPalette palette;
    palette.setColor(QPalette::Window, QColor(30, 30, 30));
    palette.setColor(QPalette::WindowText, Qt::white);
    palette.setColor(QPalette::Base, QColor(20, 20, 20));
    palette.setColor(QPalette::AlternateBase, QColor(30, 30, 30));
    palette.setColor(QPalette::ToolTipBase, Qt::white);
    palette.setColor(QPalette::ToolTipText, Qt::white);
    palette.setColor(QPalette::Text, Qt::white);
    palette.setColor(QPalette::Button, QColor(45, 45, 45));
    palette.setColor(QPalette::ButtonText, Qt::white);
    palette.setColor(QPalette::BrightText, Qt::red);
    palette.setColor(QPalette::Highlight, QColor(64, 128, 255));
    palette.setColor(QPalette::HighlightedText, Qt::black);
    app.setPalette(palette);

    app.setStyleSheet(appStyleSheet());

    MainWindow w;
    w.show();
    return app.exec();
}

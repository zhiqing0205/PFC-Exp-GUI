#include "MainWindow.h"

#include <QApplication>
#include <QIcon>
#include <QLinearGradient>
#include <QPainter>
#include <QPixmap>

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

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    app.setApplicationName("MID Nano");
#ifdef MID_NANO_VERSION
    app.setApplicationVersion(QString::fromUtf8(MID_NANO_VERSION));
#endif
    app.setWindowIcon(makeAppIcon());

    MainWindow w;
    w.show();
    return app.exec();
}

#include "PlotWidget.h"

#include <QPainter>
#include <QPainterPath>

#include <algorithm>
#include <cmath>
#include <limits>

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent) {
    setAutoFillBackground(false);
}

QSize PlotWidget::minimumSizeHint() const {
    return QSize(320, 220);
}

void PlotWidget::clear() {
    xs_.clear();
    ys_.clear();
    yLabel_.clear();
    update();
}

void PlotWidget::setData(QVector<double> xs, QVector<double> ys, QString yLabel) {
    xs_ = std::move(xs);
    ys_ = std::move(ys);
    yLabel_ = std::move(yLabel);
    update();
}

static bool isFiniteValue(double v) {
    return std::isfinite(v);
}

void PlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, true);

    const QRect r = rect();
    p.fillRect(r, QColor(14, 14, 14));

    const int left = 56;
    const int top = 18;
    const int right = 16;
    const int bottom = 34;
    const QRect plot = r.adjusted(left, top, -right, -bottom);

    p.setPen(QPen(QColor(42, 42, 42), 1));
    p.drawRoundedRect(plot.adjusted(-1, -1, 1, 1), 10, 10);

    if (xs_.isEmpty() || ys_.isEmpty()) {
        p.setPen(QColor(170, 170, 170));
        p.drawText(plot, Qt::AlignCenter, "No data");
        return;
    }

    const int n = std::min(xs_.size(), ys_.size());
    if (n <= 1) {
        p.setPen(QColor(170, 170, 170));
        p.drawText(plot, Qt::AlignCenter, "Not enough points");
        return;
    }

    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < n; ++i) {
        const double x = xs_[i];
        const double y = ys_[i];
        if (!isFiniteValue(x) || !isFiniteValue(y)) continue;
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
    }

    if (!isFiniteValue(xmin) || !isFiniteValue(xmax) || !isFiniteValue(ymin) || !isFiniteValue(ymax)) {
        p.setPen(QColor(170, 170, 170));
        p.drawText(plot, Qt::AlignCenter, "No finite values");
        return;
    }

    if (xmin == xmax) {
        xmin -= 1.0;
        xmax += 1.0;
    }
    if (ymin == ymax) {
        ymin -= 1.0;
        ymax += 1.0;
    }

    const auto mapX = [&](double x) -> double { return plot.left() + (x - xmin) / (xmax - xmin) * plot.width(); };
    const auto mapY = [&](double y) -> double { return plot.bottom() - (y - ymin) / (ymax - ymin) * plot.height(); };

    // Grid
    p.setPen(QPen(QColor(32, 32, 32), 1));
    const int gridN = 5;
    for (int i = 1; i < gridN; ++i) {
        const double t = static_cast<double>(i) / gridN;
        const int gx = plot.left() + static_cast<int>(t * plot.width());
        const int gy = plot.top() + static_cast<int>(t * plot.height());
        p.drawLine(gx, plot.top(), gx, plot.bottom());
        p.drawLine(plot.left(), gy, plot.right(), gy);
    }

    // Axes
    p.setPen(QPen(QColor(58, 58, 58), 1));
    p.drawLine(plot.bottomLeft(), plot.bottomRight());
    p.drawLine(plot.bottomLeft(), plot.topLeft());

    // Labels
    p.setPen(QColor(180, 180, 180));
    const QString xText = QString::number(xmin, 'g', 6) + " … " + QString::number(xmax, 'g', 6);
    const QString yText = QString::number(ymin, 'g', 6) + " … " + QString::number(ymax, 'g', 6);
    p.drawText(QRect(r.left(), plot.bottom() + 6, r.width(), bottom - 6), Qt::AlignCenter, xText);
    p.save();
    p.translate(12, plot.center().y());
    p.rotate(-90);
    const QString yLabel = yLabel_.isEmpty() ? "y" : yLabel_;
    p.drawText(QRect(-plot.height() / 2, -left + 10, plot.height(), left - 20), Qt::AlignCenter,
               yLabel + " (" + yText + ")");
    p.restore();

    // Curve
    QPainterPath path;
    bool started = false;
    for (int i = 0; i < n; ++i) {
        const double x = xs_[i];
        const double y = ys_[i];
        if (!isFiniteValue(x) || !isFiniteValue(y)) continue;
        const QPointF pt(mapX(x), mapY(y));
        if (!started) {
            path.moveTo(pt);
            started = true;
        } else {
            path.lineTo(pt);
        }
    }

    p.setPen(QPen(QColor(64, 128, 255), 2));
    p.drawPath(path);
}

#pragma once

#include <QVector>
#include <QWidget>

class PlotWidget final : public QWidget {
    Q_OBJECT

public:
    explicit PlotWidget(QWidget* parent = nullptr);

    void clear();
    void setData(QVector<double> xs, QVector<double> ys, QString yLabel = {});

    QSize minimumSizeHint() const override;

protected:
    void paintEvent(QPaintEvent* event) override;

private:
    QVector<double> xs_;
    QVector<double> ys_;
    QString yLabel_;
};


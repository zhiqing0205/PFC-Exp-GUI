#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QtQuickControls2/QQuickStyle>
#include <QUrl>

#include "ExperimentController.h"

int main(int argc, char* argv[]) {
  QGuiApplication app(argc, argv);

  QQuickStyle::setStyle("Material");

  ExperimentController controller;

  QQmlApplicationEngine engine;
  engine.rootContext()->setContextProperty("controller", &controller);
  engine.load(QUrl(QStringLiteral("qrc:/ExperimentGui/src/gui/qml/Main.qml")));

  if (engine.rootObjects().isEmpty()) return 1;
  return app.exec();
}

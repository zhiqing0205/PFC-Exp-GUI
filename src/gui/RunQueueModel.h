#pragma once

#include <QAbstractListModel>
#include <QDateTime>
#include <QString>
#include <QVector>

class RunQueueModel final : public QAbstractListModel {
  Q_OBJECT

public:
  enum Role {
    NameRole = Qt::UserRole + 1,
    StatusRole,
    OutDirRole,
    CommandRole,
    ExitCodeRole,
    LogTailRole,
    StartedAtRole,
    FinishedAtRole,
  };
  Q_ENUM(Role)

  struct Item {
    QString name;
    QString status;
    QString outDir;
    QString command;
    int exitCode = -1;
    QString logTail;
    QDateTime startedAt;
    QDateTime finishedAt;
  };

  explicit RunQueueModel(QObject* parent = nullptr);

  int rowCount(const QModelIndex& parent = QModelIndex()) const override;
  QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
  QHash<int, QByteArray> roleNames() const override;

  void clear();
  int addItem(const Item& item);
  const Item* itemAt(int row) const;

  void setStatus(int row, const QString& status);
  void setCommand(int row, const QString& command);
  void setOutDir(int row, const QString& outDir);
  void markStarted(int row);
  void markFinished(int row, int exitCode);
  void appendLog(int row, const QString& text);

private:
  static QString trimTail(const QString& s);

  QVector<Item> items_;
};


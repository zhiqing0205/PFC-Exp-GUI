#include "RunQueueModel.h"

#include <QVariant>

RunQueueModel::RunQueueModel(QObject* parent) : QAbstractListModel(parent) {}

int RunQueueModel::rowCount(const QModelIndex& parent) const {
  if (parent.isValid()) return 0;
  return items_.size();
}

QVariant RunQueueModel::data(const QModelIndex& index, int role) const {
  if (!index.isValid()) return {};
  const int row = index.row();
  if (row < 0 || row >= items_.size()) return {};
  const Item& item = items_.at(row);

  switch (role) {
    case NameRole:
      return item.name;
    case StatusRole:
      return item.status;
    case OutDirRole:
      return item.outDir;
    case CommandRole:
      return item.command;
    case ExitCodeRole:
      return item.exitCode;
    case LogTailRole:
      return item.logTail;
    case StartedAtRole:
      return item.startedAt.isValid() ? item.startedAt.toString(Qt::ISODate) : QString();
    case FinishedAtRole:
      return item.finishedAt.isValid() ? item.finishedAt.toString(Qt::ISODate) : QString();
    default:
      return {};
  }
}

QHash<int, QByteArray> RunQueueModel::roleNames() const {
  return {
      {NameRole, "name"},
      {StatusRole, "status"},
      {OutDirRole, "outDir"},
      {CommandRole, "command"},
      {ExitCodeRole, "exitCode"},
      {LogTailRole, "logTail"},
      {StartedAtRole, "startedAt"},
      {FinishedAtRole, "finishedAt"},
  };
}

void RunQueueModel::clear() {
  beginResetModel();
  items_.clear();
  endResetModel();
}

int RunQueueModel::addItem(const Item& item) {
  const int row = items_.size();
  beginInsertRows(QModelIndex(), row, row);
  items_.push_back(item);
  endInsertRows();
  return row;
}

const RunQueueModel::Item* RunQueueModel::itemAt(int row) const {
  if (row < 0 || row >= items_.size()) return nullptr;
  return &items_.at(row);
}

void RunQueueModel::setStatus(int row, const QString& status) {
  if (row < 0 || row >= items_.size()) return;
  if (items_[row].status == status) return;
  items_[row].status = status;
  const QModelIndex idx = index(row);
  emit dataChanged(idx, idx, {StatusRole});
}

void RunQueueModel::setCommand(int row, const QString& command) {
  if (row < 0 || row >= items_.size()) return;
  if (items_[row].command == command) return;
  items_[row].command = command;
  const QModelIndex idx = index(row);
  emit dataChanged(idx, idx, {CommandRole});
}

void RunQueueModel::setOutDir(int row, const QString& outDir) {
  if (row < 0 || row >= items_.size()) return;
  if (items_[row].outDir == outDir) return;
  items_[row].outDir = outDir;
  const QModelIndex idx = index(row);
  emit dataChanged(idx, idx, {OutDirRole});
}

void RunQueueModel::markStarted(int row) {
  if (row < 0 || row >= items_.size()) return;
  items_[row].startedAt = QDateTime::currentDateTime();
  const QModelIndex idx = index(row);
  emit dataChanged(idx, idx, {StartedAtRole});
}

void RunQueueModel::markFinished(int row, int exitCode) {
  if (row < 0 || row >= items_.size()) return;
  items_[row].exitCode = exitCode;
  items_[row].finishedAt = QDateTime::currentDateTime();
  const QModelIndex idx = index(row);
  emit dataChanged(idx, idx, {ExitCodeRole, FinishedAtRole});
}

void RunQueueModel::appendLog(int row, const QString& text) {
  if (row < 0 || row >= items_.size()) return;
  items_[row].logTail = trimTail(items_[row].logTail + text);
  const QModelIndex idx = index(row);
  emit dataChanged(idx, idx, {LogTailRole});
}

QString RunQueueModel::trimTail(const QString& s) {
  static constexpr int kMaxChars = 12000;
  if (s.size() <= kMaxChars) return s;
  return s.right(kMaxChars);
}


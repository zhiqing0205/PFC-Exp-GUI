import QtQuick
import QtQuick.Controls
import QtQuick.Dialogs
import QtQuick.Layouts

ApplicationWindow {
    id: root
    width: 1180
    height: 760
    visible: true
    title: "实验参数配置与批量运行"

    function numOr(text, fallback) {
        const v = Number(text)
        return isNaN(v) ? fallback : v
    }

    function intOr(text, fallback) {
        const v = Number(text)
        if (isNaN(v)) return fallback
        return Math.max(1, Math.floor(v))
    }

    function countRange(startText, endText, stepText) {
        const start = Number(startText)
        const end = Number(endText)
        const step = Number(stepText)
        if (isNaN(start) || isNaN(end) || isNaN(step)) return 0
        if (step <= 0) return 0
        if (end < start) return 0
        return Math.floor((end - start) / step + 1.0000000001)
    }

    function previewRunCount() {
        const cU0 = u0Sweep.checked ? countRange(u0Start.text, u0End.text, u0Step.text) : 1
        const cCon0 = con0Sweep.checked ? countRange(con0Start.text, con0End.text, con0Step.text) : 1
        const cSteps = stepsSweep.checked ? countRange(stepsStart.value, stepsEnd.value, stepsStep.value) : 1
        return cU0 * cCon0 * cSteps
    }

    function buildConfig() {
        return {
            simPath: simPathField.text,
            useMpi: useMpiCheck.checked,
            mpiLauncher: mpiLauncherField.text,
            mpiRanks: mpiRanksSpin.value,
            outBase: outBaseField.text,
            runPrefix: runPrefixField.text,
            params: {
                u0: numOr(u0Single.text, 0.05),
                con0: numOr(con0Single.text, 0.2),
                sig: numOr(sigField.text, 0.05),
                dt: numOr(dtField.text, 0.05),
                dx: numOr(dxField.text, 0.125),
                steps: intOr(stepsSingle.value, 5002),
                mod: intOr(modField.value, 25),
                grainx: intOr(grainxField.value, 80),
                grainy: intOr(grainyField.value, 80),
                grainz: intOr(grainzField.value, 1),
                axx: Math.max(0, numOr(axxField.text, 0.0)),
            },
            sweeps: {
                u0: {
                    enabled: u0Sweep.checked,
                    start: numOr(u0Start.text, numOr(u0Single.text, 0.05)),
                    end: numOr(u0End.text, numOr(u0Single.text, 0.05)),
                    step: numOr(u0Step.text, 1.0),
                },
                con0: {
                    enabled: con0Sweep.checked,
                    start: numOr(con0Start.text, numOr(con0Single.text, 0.2)),
                    end: numOr(con0End.text, numOr(con0Single.text, 0.2)),
                    step: numOr(con0Step.text, 1.0),
                },
                steps: {
                    enabled: stepsSweep.checked,
                    start: stepsStart.value,
                    end: stepsEnd.value,
                    step: stepsStep.value,
                },
            }
        }
    }

    header: ToolBar {
        RowLayout {
            anchors.fill: parent
            spacing: 12
            Label {
                text: root.title
                font.pixelSize: 18
                Layout.fillWidth: true
            }
            Label {
                text: controller.running ? ("运行中：" + controller.finishedRuns + "/" + controller.totalRuns) : ("已完成：" + controller.finishedRuns + "/" + controller.totalRuns)
                opacity: 0.85
            }
        }
    }

    SplitView {
        anchors.fill: parent

        ScrollView {
            SplitView.preferredWidth: 560
            contentWidth: availableWidth

            ColumnLayout {
                width: parent.width
                spacing: 12
                padding: 12

                GroupBox {
                    title: "程序与运行方式"
                    Layout.fillWidth: true
                    GridLayout {
                        columns: 3
                        columnSpacing: 8
                        rowSpacing: 8
                        anchors.fill: parent

                        Label { text: "实验程序 (misfit_sim)"; Layout.alignment: Qt.AlignVCenter }
                        TextField { id: simPathField; placeholderText: "留空则默认同目录"; Layout.fillWidth: true }
                        Button {
                            text: "选择…"
                            onClicked: simDialog.open()
                        }

                        Label { text: "使用 MPI"; Layout.alignment: Qt.AlignVCenter }
                        CheckBox { id: useMpiCheck; text: ""; checked: false }
                        Item { }

                        Label { text: "MPI 启动器"; Layout.alignment: Qt.AlignVCenter; enabled: useMpiCheck.checked }
                        TextField { id: mpiLauncherField; text: "mpirun"; enabled: useMpiCheck.checked; Layout.fillWidth: true }
                        Item { }

                        Label { text: "MPI 进程数"; Layout.alignment: Qt.AlignVCenter; enabled: useMpiCheck.checked }
                        SpinBox { id: mpiRanksSpin; from: 1; to: 1024; value: 1; enabled: useMpiCheck.checked; Layout.fillWidth: true }
                        Item { }
                    }
                }

                GroupBox {
                    title: "输出目录"
                    Layout.fillWidth: true
                    GridLayout {
                        columns: 3
                        columnSpacing: 8
                        rowSpacing: 8
                        anchors.fill: parent

                        Label { text: "父目录"; Layout.alignment: Qt.AlignVCenter }
                        TextField { id: outBaseField; text: "runs"; Layout.fillWidth: true }
                        Button { text: "选择…"; onClicked: outDialog.open() }

                        Label { text: "会话前缀"; Layout.alignment: Qt.AlignVCenter }
                        TextField { id: runPrefixField; text: "session"; Layout.fillWidth: true }
                        Item { }

                        Label { text: "预估组合数"; Layout.alignment: Qt.AlignVCenter }
                        Label { text: previewRunCount(); font.bold: true; Layout.fillWidth: true }
                        Item { }
                    }
                }

                GroupBox {
                    title: "关键参数（单次/扫描）"
                    Layout.fillWidth: true

                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 10

                        GridLayout {
                            columns: 6
                            columnSpacing: 8
                            rowSpacing: 6
                            Layout.fillWidth: true

                            Label { text: ""; }
                            Label { text: "单值"; }
                            Label { text: "扫描"; }
                            Label { text: "起始"; }
                            Label { text: "结束"; }
                            Label { text: "步长"; }

                            Label { text: "密度 u0"; Layout.alignment: Qt.AlignVCenter }
                            TextField { id: u0Single; text: "0.05"; inputMethodHints: Qt.ImhFormattedNumbersOnly }
                            CheckBox { id: u0Sweep; text: "" }
                            TextField { id: u0Start; text: "0.05"; enabled: u0Sweep.checked; inputMethodHints: Qt.ImhFormattedNumbersOnly }
                            TextField { id: u0End; text: "0.05"; enabled: u0Sweep.checked; inputMethodHints: Qt.ImhFormattedNumbersOnly }
                            TextField { id: u0Step; text: "0.01"; enabled: u0Sweep.checked; inputMethodHints: Qt.ImhFormattedNumbersOnly }

                            Label { text: "浓度 con0"; Layout.alignment: Qt.AlignVCenter }
                            TextField { id: con0Single; text: "0.2"; inputMethodHints: Qt.ImhFormattedNumbersOnly }
                            CheckBox { id: con0Sweep; text: "" }
                            TextField { id: con0Start; text: "0.2"; enabled: con0Sweep.checked; inputMethodHints: Qt.ImhFormattedNumbersOnly }
                            TextField { id: con0End; text: "0.2"; enabled: con0Sweep.checked; inputMethodHints: Qt.ImhFormattedNumbersOnly }
                            TextField { id: con0Step; text: "0.05"; enabled: con0Sweep.checked; inputMethodHints: Qt.ImhFormattedNumbersOnly }

                            Label { text: "迭代步数 steps"; Layout.alignment: Qt.AlignVCenter }
                            SpinBox { id: stepsSingle; from: 1; to: 2000000; value: 5002 }
                            CheckBox { id: stepsSweep; text: "" }
                            SpinBox { id: stepsStart; from: 1; to: 2000000; value: 5002; enabled: stepsSweep.checked }
                            SpinBox { id: stepsEnd; from: 1; to: 2000000; value: 5002; enabled: stepsSweep.checked }
                            SpinBox { id: stepsStep; from: 1; to: 2000000; value: 100; enabled: stepsSweep.checked }
                        }
                    }
                }

                GroupBox {
                    title: "其它参数（单值）"
                    Layout.fillWidth: true
                    GridLayout {
                        columns: 4
                        columnSpacing: 8
                        rowSpacing: 8
                        anchors.fill: parent

                        Label { text: "sig"; Layout.alignment: Qt.AlignVCenter }
                        TextField { id: sigField; text: "0.05"; inputMethodHints: Qt.ImhFormattedNumbersOnly; Layout.fillWidth: true }
                        Label { text: "dt"; Layout.alignment: Qt.AlignVCenter }
                        TextField { id: dtField; text: "0.05"; inputMethodHints: Qt.ImhFormattedNumbersOnly; Layout.fillWidth: true }

                        Label { text: "dx"; Layout.alignment: Qt.AlignVCenter }
                        TextField { id: dxField; text: "0.125"; inputMethodHints: Qt.ImhFormattedNumbersOnly; Layout.fillWidth: true }
                        Label { text: "mod"; Layout.alignment: Qt.AlignVCenter }
                        SpinBox { id: modField; from: 1; to: 1000000; value: 25; Layout.fillWidth: true }

                        Label { text: "grainx"; Layout.alignment: Qt.AlignVCenter }
                        SpinBox { id: grainxField; from: 1; to: 100000; value: 80; Layout.fillWidth: true }
                        Label { text: "grainy"; Layout.alignment: Qt.AlignVCenter }
                        SpinBox { id: grainyField; from: 1; to: 100000; value: 80; Layout.fillWidth: true }

                        Label { text: "grainz"; Layout.alignment: Qt.AlignVCenter }
                        SpinBox { id: grainzField; from: 1; to: 100000; value: 1; Layout.fillWidth: true }
                        Label { text: "噪声 axx"; Layout.alignment: Qt.AlignVCenter }
                        TextField { id: axxField; text: "0"; inputMethodHints: Qt.ImhFormattedNumbersOnly; Layout.fillWidth: true }
                    }
                }

                RowLayout {
                    Layout.fillWidth: true
                    spacing: 8
                    Button {
                        text: "开始"
                        enabled: !controller.running
                        onClicked: controller.startRuns(buildConfig())
                    }
                    Button {
                        text: "停止"
                        enabled: controller.running
                        onClicked: controller.stop()
                    }
                    Button {
                        text: "清空"
                        enabled: !controller.running
                        onClicked: controller.clear()
                    }
                    Item { Layout.fillWidth: true }
                }

                Label {
                    Layout.fillWidth: true
                    visible: controller.lastError.length > 0
                    color: "#C62828"
                    text: controller.lastError
                    wrapMode: Text.Wrap
                }
            }
        }

        Rectangle {
            color: Qt.rgba(0, 0, 0, 0)
            SplitView.fillWidth: true

            ColumnLayout {
                anchors.fill: parent
                spacing: 8
                padding: 12

                GroupBox {
                    title: "任务队列"
                    Layout.fillWidth: true
                    Layout.preferredHeight: 300
                    ListView {
                        id: runList
                        anchors.fill: parent
                        clip: true
                        model: controller.runModel
                        delegate: ItemDelegate {
                            width: ListView.view.width
                            text: name + "  [" + status + "]"
                            onClicked: controller.selectedRun = index
                        }
                    }
                }

                GroupBox {
                    title: "日志（尾部）"
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    TextArea {
                        anchors.fill: parent
                        readOnly: true
                        wrapMode: Text.Wrap
                        text: controller.selectedLog
                        font.family: "monospace"
                    }
                }
            }
        }
    }

    FileDialog {
        id: simDialog
        title: "选择实验程序 (misfit_sim)"
        onAccepted: simPathField.text = selectedFile.toLocalFile()
    }

    FolderDialog {
        id: outDialog
        title: "选择输出父目录"
        onAccepted: outBaseField.text = selectedFolder.toLocalFile()
    }
}


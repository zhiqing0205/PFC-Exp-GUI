#include <windows.h>

#include <shellapi.h>

#include <string>
#include <vector>

static std::wstring exeDir() {
    std::wstring path;
    path.resize(MAX_PATH);
    DWORD len = GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
    if (len == 0) return L".";
    path.resize(len);

    const size_t pos = path.find_last_of(L"\\/");
    if (pos == std::wstring::npos) return L".";
    return path.substr(0, pos);
}

static std::wstring quoteArg(const std::wstring& arg) {
    if (arg.empty()) return L"\"\"";

    bool needsQuotes = false;
    for (wchar_t c : arg) {
        if (c == L' ' || c == L'\t' || c == L'\n' || c == L'\v' || c == L'\"') {
            needsQuotes = true;
            break;
        }
    }
    if (!needsQuotes) return arg;

    std::wstring out;
    out.push_back(L'"');
    for (wchar_t c : arg) {
        if (c == L'"') out.push_back(L'\\');
        out.push_back(c);
    }
    out.push_back(L'"');
    return out;
}

static void showError(const std::wstring& message) {
    MessageBoxW(nullptr, message.c_str(), L"MID Nano", MB_ICONERROR | MB_OK);
}

int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR, int) {
    const std::wstring dir = exeDir();
    const std::wstring appDir = dir + L"\\app";
    const std::wstring targetExe = appDir + L"\\MID Nano.exe";

    if (GetFileAttributesW(targetExe.c_str()) == INVALID_FILE_ATTRIBUTES) {
        showError(L"Cannot find the packaged app executable:\n" + targetExe + L"\n\nRe-download the ZIP release and keep the 'app' folder next to MID Nano.exe.");
        return 1;
    }

    int argc = 0;
    LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);

    std::wstring cmdLine = quoteArg(targetExe);
    if (argv) {
        for (int i = 1; i < argc; ++i) {
            cmdLine.push_back(L' ');
            cmdLine += quoteArg(argv[i] ? std::wstring(argv[i]) : std::wstring());
        }
        LocalFree(argv);
    }

    STARTUPINFOW si{};
    si.cb = sizeof(si);
    PROCESS_INFORMATION pi{};

    std::vector<wchar_t> mutableCmd(cmdLine.begin(), cmdLine.end());
    mutableCmd.push_back(L'\0');

    const BOOL ok = CreateProcessW(
        targetExe.c_str(),
        mutableCmd.data(),
        nullptr,
        nullptr,
        FALSE,
        0,
        nullptr,
        appDir.c_str(),
        &si,
        &pi);
    if (!ok) {
        showError(L"Failed to launch the packaged app.\n\nError code: " + std::to_wstring(GetLastError()));
        return 2;
    }

    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    return 0;
}


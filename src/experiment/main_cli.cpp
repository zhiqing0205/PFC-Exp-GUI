// main_cli.cpp - Unified CLI entry point for PFC simulations.
// Dispatches to the appropriate model based on --model argument.
#include <cstring>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

// Model entry points (defined in separate .cpp files).
int run_misfit(int argc, char** argv);
int run_cvd(int argc, char** argv);
int run_elastic(int argc, char** argv);

static void PrintGlobalUsage(const char* argv0) {
    std::cerr
        << "Usage: " << (argv0 ? argv0 : "pfc-exp-cli") << " [--model <name>] [model-options...]\n\n"
        << "Models:\n"
        << "  misfit    Phase-field crystal with misfit strain (default)\n"
        << "  cvd       Controlled vapor deposition / two-phase growth\n"
        << "  elastic   Post-processing: elastic strain & crystallographic order\n\n"
        << "Use --model <name> --help for model-specific options.\n";
}

int main(int argc, char** argv) {
    std::string model_name = "misfit";  // default

    // Scan for --model and extract it, building a new argv without --model.
    std::vector<char*> new_argv;
    new_argv.push_back(argv[0]);

    for (int i = 1; i < argc; ++i) {
        if (!argv[i]) continue;
        std::string_view arg(argv[i]);

        if (arg == "--model" || arg == "-m") {
            if (i + 1 < argc) {
                model_name = argv[++i];
            } else {
                std::cerr << "Missing value after --model\n";
                PrintGlobalUsage(argv[0]);
                return 2;
            }
            continue;
        }

        // Handle --model=name
        if (arg.substr(0, 8) == "--model=") {
            model_name = std::string(arg.substr(8));
            continue;
        }

        // Global help (before model dispatch)
        if (arg == "-h" || arg == "--help") {
            // If no --model was specified yet, show global help
            // Otherwise, pass --help to the model
            if (model_name == "misfit" && i == 1) {
                // Check if there's a --model later
                bool has_model = false;
                for (int j = i + 1; j < argc; ++j) {
                    if (argv[j] && (std::string_view(argv[j]) == "--model" || std::string_view(argv[j]).substr(0, 8) == "--model=")) {
                        has_model = true;
                        break;
                    }
                }
                if (!has_model) {
                    PrintGlobalUsage(argv[0]);
                    return 0;
                }
            }
        }

        new_argv.push_back(argv[i]);
    }

    int new_argc = static_cast<int>(new_argv.size());

    if (model_name == "misfit") {
        return run_misfit(new_argc, new_argv.data());
    } else if (model_name == "cvd") {
        return run_cvd(new_argc, new_argv.data());
    } else if (model_name == "elastic") {
        return run_elastic(new_argc, new_argv.data());
    } else {
        std::cerr << "Unknown model: " << model_name << "\n";
        PrintGlobalUsage(argv[0]);
        return 2;
    }
}

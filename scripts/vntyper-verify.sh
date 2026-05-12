#!/usr/bin/env bash
# Verify the Rust kestrel binary is usable through VNtyper's Kestrel path.
#
# This installs the Rust-built binary at VNtyper's binary dependency location:
#   kestrel/integration/VNtyper/vntyper/dependencies/kestrel/kestrel
# and creates a small compatibility kestrel.jar in the same directory. VNtyper
# still constructs `java -jar .../kestrel.jar ...`; the compatibility jar execs
# the sibling Rust binary, so the unchanged VNtyper command shape reaches Rust.
#
# The VNtyper unit suite does not run the full BAM/FASTQ pipeline or require the
# large Zenodo integration data set. It validates VNtyper's import/configuration
# surface while this script separately probes the Java-shaped VNtyper Kestrel
# invocation and confirms it reaches the Rust CLI.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

KESTREL_BIN="${1:-${KESTREL_BIN:-target/release/kestrel}}"
VNTYPER="$REPO_ROOT/kestrel/integration/VNtyper"
KESTREL_DEPS="$VNTYPER/vntyper/dependencies/kestrel"
VENV="${VENV:-$REPO_ROOT/.venv-vntyper}"

if [[ ! -x "$KESTREL_BIN" ]]; then
    echo "FAIL: Rust kestrel binary is missing or not executable: $KESTREL_BIN" >&2
    echo "Build it with: cargo build --release -p kestrel --bin kestrel" >&2
    exit 1
fi

if [[ ! -d "$VNTYPER" ]]; then
    echo "FAIL: VNtyper submodule missing at $VNTYPER" >&2
    echo "Run: git submodule update --init --recursive kestrel" >&2
    exit 1
fi

echo ">> Installing Rust kestrel binary into $KESTREL_DEPS"
mkdir -p "$KESTREL_DEPS"
cp "$KESTREL_BIN" "$KESTREL_DEPS/kestrel"
chmod +x "$KESTREL_DEPS/kestrel"

if ! command -v javac >/dev/null || ! command -v jar >/dev/null; then
    echo "FAIL: javac and jar are required to build the VNtyper compatibility launcher" >&2
    exit 1
fi

echo ">> Building compatibility kestrel.jar that delegates to the Rust binary"
WRAPPER_DIR="$(mktemp -d)"
trap 'rm -rf "$WRAPPER_DIR"' EXIT
cat > "$WRAPPER_DIR/RustKestrelLauncher.java" <<'JAVA'
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public final class RustKestrelLauncher {
    public static void main(String[] args) throws Exception {
        Path jarPath = Path.of(RustKestrelLauncher.class
            .getProtectionDomain()
            .getCodeSource()
            .getLocation()
            .toURI());
        Path rustBinary = jarPath.getParent().resolve("kestrel");

        List<String> command = new ArrayList<>();
        command.add(rustBinary.toString());
        for (String arg : args) {
            command.add(arg);
        }

        Process process = new ProcessBuilder(command)
            .inheritIO()
            .start();
        System.exit(process.waitFor());
    }
}
JAVA
cat > "$WRAPPER_DIR/MANIFEST.MF" <<'EOF'
Manifest-Version: 1.0
Main-Class: RustKestrelLauncher

EOF
javac -d "$WRAPPER_DIR/classes" "$WRAPPER_DIR/RustKestrelLauncher.java"
jar cfm "$KESTREL_DEPS/kestrel.jar" "$WRAPPER_DIR/MANIFEST.MF" -C "$WRAPPER_DIR/classes" .

if [[ ! -x "$VENV/bin/python" ]]; then
    echo ">> Creating Python virtualenv at $VENV"
    python3 -m venv "$VENV"
fi

echo ">> Installing VNtyper test dependencies"
"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install -e "$VNTYPER[dev]"

echo ">> Running VNtyper unit tests"
(
    cd "$VNTYPER"
    "$VENV/bin/python" -m pytest tests/unit/ -q --no-header
)

echo ""
echo ">> Probing VNtyper's Java-shaped Kestrel invocation"
if ! java -jar "$KESTREL_DEPS/kestrel.jar" -h | grep -q "Rust port of the Kestrel variant caller"; then
    echo "FAIL: java -jar $KESTREL_DEPS/kestrel.jar did not reach the Rust kestrel CLI" >&2
    exit 1
fi

echo ""
echo "PASS: VNtyper unit tests passed and java -jar kestrel.jar reaches Rust kestrel"

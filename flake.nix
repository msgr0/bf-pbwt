{
  description = "MacOS C development environment with GCC and HTSlib";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        
        # On macOS, pkgs.stdenv is Clang-based. 
        # We explicitly use gccStdenv to ensure 'cc' points to GCC.

      in
      {
        devShells.default = pkgs.mkShell.override {} {
          buildInputs = with pkgs; [
          bcftools
            clang-tools
            clang
            llvmPackages.openmp
            libcxx
            libunwind
            cmake
            codespell
            conan
            ninja
            doxygen
            htslib
            pkg-config
            zlib
            bzip2
            xz
            curl
            bear
          ] ++ (if system == "aarch64-darwin" then [ ] else [ gdb ]);  

          shellHook = ''
            echo "--- macOS GCC Dev Environment ---"
            echo "Compiler: $(cc --version | head -n 1)"
            echo "HTSlib version: $(pkg-config --modversion htslib)"
            export CFLAGS="-I${pkgs.llvmPackages.openmp}/include $CFLAGS"
            export LDFLAGS="-L${pkgs.llvmPackages.openmp}/lib -lomp $LDFLAGS"
          '';
        };
      });
}

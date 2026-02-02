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
        stdenv = pkgs.gccStdenv;
      in
      {
        devShells.default = pkgs.mkShell.override { inherit stdenv; } {
          buildInputs = with pkgs; [
            gcc
            
            # Bioinformatics
            htslib
           
            pkg-config
            
            # Common HTSlib dependencies (often needed for linking)
            zlib
            bzip2
            xz
            curl
            bear
          ];

          shellHook = ''
            echo "--- macOS GCC Dev Environment ---"
            echo "Compiler: $(cc --version | head -n 1)"
            echo "HTSlib version: $(pkg-config --modversion htslib)"
          '';
        };
      });
}

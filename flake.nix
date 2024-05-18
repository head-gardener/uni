{
  description = "Description for the project";

  inputs = {
    flake-parts.url = "github:hercules-ci/flake-parts";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      imports = [

      ];
      systems = [ "x86_64-linux" ];
      perSystem = { config, self', inputs', pkgs, system, ... }: {

        devShells.default = pkgs.mkShell {
          name = "edu";
          packages = with pkgs; [
            just
            julia
          ];
          shellHook = ''
          echo Have fun!
          '';
        };
      };
      flake = {

      };
    };
}

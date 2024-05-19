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
      perSystem = { config, self', inputs', pkgs, lib, system, ... }: 
      let
        # https://s3.backyard-hg.xyz/stuff/CiscoPacketTracer_821_Ubuntu_64bit.deb
        cpt821 = (pkgs.ciscoPacketTracer8.overrideAttrs (_: _: { version = "8.2.1"; }));
        pt = pkgs.writeScriptBin "pt" ''
          unshare -rn ${cpt821}/bin/packettracer8
        '';
      in {
        _module.args.pkgs = import inputs.nixpkgs {
          inherit system;
          config.allowUnfreePredicate = pkg: builtins.elem (lib.getName pkg) [
            "ciscoPacketTracer8"
          ];
        };

        devShells.default = pkgs.mkShell {
          name = "edu";
          packages = with pkgs; [
            pt
            just
            julia
          ];
          shellHook = ''
            echo Have fun!
          '';
        };
      };
      flake = { };
    };
}

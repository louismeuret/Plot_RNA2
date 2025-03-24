import { EmptyLoci } from "../../mol-model/loci";
import { StructureSelection } from "../../mol-model/structure";
import { BuiltInTrajectoryFormat } from "../../mol-plugin-state/formats/trajectory";
import { createStructureRepresentationParams } from "../../mol-plugin-state/helpers/structure-representation-params";
import { setStructureOverpaint } from "../../mol-plugin-state/helpers/structure-overpaint";
import { createPluginUI } from "../../mol-plugin-ui";
import { PluginUIContext } from "../../mol-plugin-ui/context";
import { renderReact18 } from "../../mol-plugin-ui/react18";
import { DefaultPluginUISpec } from "../../mol-plugin-ui/spec";
import { PluginCommands } from "../../mol-plugin/commands";
import { StateTransforms } from "../../mol-plugin-state/transforms";
import { StateElements } from "./helpers";
import { StateObject } from "../../mol-state";
import { Script } from "../../mol-script/script";
import { StateObjectSelector } from "../../mol-state";
import { Color } from "../../mol-util/color";
import { StripedResidues } from "./coloring";
import { CustomColorThemeProvider } from "./custom-theme";
import "./index.html";
import "../../mol-plugin-ui/skin/light.scss";
import { PluginStateObject } from "../../mol-plugin-state/objects";
import { TrajectoryFromModelAndCoordinates } from "../../mol-plugin-state/transforms/model";
import { BuiltInTopologyFormat } from "../../mol-plugin-state/formats/topology";
import { BuiltInCoordinatesFormat } from "../../mol-plugin-state/formats/coordinates";
import { PresetTrajectoryHierarchy } from "../../mol-plugin-state/builder/structure/hierarchy-preset";

export interface LoadTrajectoryParams {
  model:
    | {
        kind: "model-url";
        url: string;
        format?: BuiltInTrajectoryFormat /* mmcif */;
        isBinary?: boolean;
      }
    | {
        kind: "model-data";
        data: string | number[] | ArrayBuffer | Uint8Array;
        format?: BuiltInTrajectoryFormat /* mmcif */;
      }
    | {
        kind: "topology-url";
        url: string;
        format: BuiltInTopologyFormat;
        isBinary?: boolean;
      }
    | {
        kind: "topology-data";
        data: string | number[] | ArrayBuffer | Uint8Array;
        format: BuiltInTopologyFormat;
      };
  modelLabel?: string;
  coordinates:
    | {
        kind: "coordinates-url";
        url: string;
        format: BuiltInCoordinatesFormat;
        isBinary?: boolean;
      }
    | {
        kind: "coordinates-data";
        data: string | number[] | ArrayBuffer | Uint8Array;
        format: BuiltInCoordinatesFormat;
      };
  coordinatesLabel?: string;
  preset?: keyof PresetTrajectoryHierarchy;
}

class BasicWrapper {
  plugin: PluginUIContext;

  async init(target: string | HTMLElement) {
    this.plugin = await createPluginUI({
      target:
        typeof target === "string" ? document.getElementById(target)! : target,
      render: renderReact18,
      spec: {
        ...DefaultPluginUISpec(),
        layout: {
          initial: {
            isExpanded: false,
            showControls: false,
          },
        },
        components: {
          remoteState: "none",
        },
      },
    });

    this.plugin.representation.structure.themes.colorThemeRegistry.add(
      StripedResidues.colorThemeProvider!,
    );
    this.plugin.representation.structure.themes.colorThemeRegistry.add(
      CustomColorThemeProvider,
    );
    this.plugin.managers.lociLabels.addProvider(StripedResidues.labelProvider!);
    this.plugin.customModelProperties.register(
      StripedResidues.propertyProvider,
      true,
    );

    this.plugin.managers.dragAndDrop.addHandler("custom-wrapper", (files) => {
      if (files.some((f) => f.name.toLowerCase().endsWith(".testext"))) {
        console.log(".testext File dropped");
        return true;
      }
      return false;
    });
  }

  async loadTrajectory(params: LoadTrajectoryParams) {
    const plugin = this.plugin;

    let model: StateObjectSelector;

    if (
      params.model.kind === "model-data" ||
      params.model.kind === "model-url"
    ) {
      const data =
        params.model.kind === "model-data"
          ? await plugin.builders.data.rawData({
              data: params.model.data,
              label: params.modelLabel,
            })
          : await plugin.builders.data.download({
              url: params.model.url,
              isBinary: params.model.isBinary,
              label: params.modelLabel,
            });

      const trajectory = await plugin.builders.structure.parseTrajectory(
        data,
        params.model.format ?? "mmcif",
      );
      model = await plugin.builders.structure.createModel(trajectory);
    } else {
      const data =
        params.model.kind === "topology-data"
          ? await plugin.builders.data.rawData({
              data: params.model.data,
              label: params.modelLabel,
            })
          : await plugin.builders.data.download({
              url: params.model.url,
              isBinary: params.model.isBinary,
              label: params.modelLabel,
            });

      const provider = plugin.dataFormats.get(params.model.format);
      model = await provider!.parse(plugin, data);
    }

    const data =
      params.coordinates.kind === "coordinates-data"
        ? await plugin.builders.data.rawData({
            data: params.coordinates.data,
            label: params.coordinatesLabel,
          })
        : await plugin.builders.data.download({
            url: params.coordinates.url,
            isBinary: params.coordinates.isBinary,
            label: params.coordinatesLabel,
          });

    const provider = plugin.dataFormats.get(params.coordinates.format);
    const coords = await provider!.parse(plugin, data);

    const trajectory = await plugin
      .build()
      .toRoot()
      .apply(
        TrajectoryFromModelAndCoordinates,
        {
          modelRef: model.ref,
          coordinatesRef: coords.ref,
        },
        { dependsOn: [model.ref, coords.ref] },
      )
      .commit();

    const preset = await plugin.builders.structure.hierarchy.applyPreset(
      trajectory,
      params.preset ?? "default",
    );

    return { model, coords, preset };
  }

  setBackground(color: number) {
    PluginCommands.Canvas3D.SetSettings(this.plugin, {
      settings: (props) => {
        props.renderer.backgroundColor = Color(color);
      },
    });
  }

  toggleSpin() {
    if (!this.plugin.canvas3d) return;

    const trackball = this.plugin.canvas3d.props.trackball;
    PluginCommands.Canvas3D.SetSettings(this.plugin, {
      settings: {
        trackball: {
          ...trackball,
          animate:
            trackball.animate.name === "spin"
              ? { name: "off", params: {} }
              : { name: "spin", params: { speed: 1 } },
        },
      },
    });
    if (this.plugin.canvas3d.props.trackball.animate.name !== "spin") {
      PluginCommands.Camera.Reset(this.plugin, {});
    }
  }

  highlightResidue(seqId: number, color = 0xff0000) {
    const data =
      this.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj
        ?.data;
    if (!data) return;
    const sel = Script.getStructureSelection(
      (Q) =>
        Q.struct.generator.atomGroups({
          "residue-test": Q.core.rel.eq([
            Q.struct.atomProperty.macromolecular.label_seq_id(),
            seqId,
          ]),
          "group-by": Q.struct.atomProperty.macromolecular.residueKey(),
        }),
      data,
    );
    const loci = StructureSelection.toLociWithSourceUnits(sel);
    this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
    //this.plugin.managers.structure.component.updateRepresentationsTheme(data.models,      { color },    );
    this.plugin.managers.camera.focusLoci(loci);
  }

  highlightResidues(seqIds: number[], color = 0xff0000) {
    seqIds.forEach((id: number) => this.highlightResidue(id, color));
  }

  applyRepresentation(type2 = "cartoon") {
    this.plugin.dataTransaction(async () => {
      for (const s of this.plugin.managers.structure.hierarchy.current
        .structures) {
        await this.plugin.managers.structure.component.updateRepresentationsTheme(
          s.components,
          { color: "default" },
        );
      }
    });
  }

  get state() {
    return this.plugin.state.data;
  }

  private getObj<T extends StateObject>(ref: string): T["data"] {
    const state = this.state;
    const cell = state.select(ref)[0];
    if (!cell || !cell.obj) return void 0;
    return (cell.obj as T).data;
  }

  toggleWater() {
    const structure = this.getObj<PluginStateObject.Molecule.Structure>(
      StateElements.Assembly,
    );
    if (!structure) return;

    const waterElement = this.getObj<PluginStateObject.Molecule.Structure>(
      StateElements.Water,
    );
    const update = this.state.build();

    if (waterElement) {
      // Water is currently visible, so remove it
      update.to(StateElements.Water).delete(StateElements.WaterVisual);
    } else {
      // Water is hidden, so add it
      update.to(StateElements.Water).applyOrUpdate(
        StateElements.WaterVisual,
        StateTransforms.Representation.StructureRepresentation3D,
        createStructureRepresentationParams(this.plugin, structure, {
          type: "ball-and-stick",
          typeParams: { alpha: 0.51 },
          color: "element-symbol",
        }),
      );
    }

    return this.state.updateTree(update);
  }

  coloring = {
    applyStripes: async () => {
      this.plugin.dataTransaction(async () => {
        for (const s of this.plugin.managers.structure.hierarchy.current
          .structures) {
          await this.plugin.managers.structure.component.updateRepresentationsTheme(
            s.components,
            { color: StripedResidues.propertyProvider.descriptor.name as any },
          );
        }
      });
    },
    applyCustomTheme: async () => {
      this.plugin.dataTransaction(async () => {
        for (const s of this.plugin.managers.structure.hierarchy.current
          .structures) {
          await this.plugin.managers.structure.component.updateRepresentationsTheme(
            s.components,
            { color: CustomColorThemeProvider.name as any },
          );
        }
      });
    },
    applyDefault: async () => {
      this.plugin.dataTransaction(async () => {
        for (const s of this.plugin.managers.structure.hierarchy.current
          .structures) {
          await this.plugin.managers.structure.component.updateRepresentationsTheme(
            s.components,
            { color: "default" },
          );
        }
      });
    },
  };

  interactivity = {
    changeFrame: async (idx: number) => {
      if (!this.plugin) {
        console.error("Viewer not initialized.");
        return;
      }

      const plugin = this.plugin;
      const state = plugin.state.data;
      const models = plugin.managers.structure.hierarchy.current.models;
      if (!models.length) {
        console.error("No models found in the hierarchy.");
        return;
      }

      let model = models[1].cell.transform.ref;
      let update = state.build();
      update.to(model).update({ modelIndex: idx });
      update.commit();

      try {
        state.updateTree(update);
        console.log(`Frame updated to ${idx}`);
      } catch (err) {
        console.error("Error updating frame:", err);
      }
    },

    highlightOn: (seqId: number, color = 0xff0000) => {
      const data =
        this.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj
          ?.data;
      if (!data) return;

      const sel = Script.getStructureSelection(
        (Q) =>
          Q.struct.generator.atomGroups({
            "residue-test": Q.core.rel.eq([
              Q.struct.atomProperty.macromolecular.label_seq_id(),
              seqId,
            ]),
            "group-by": Q.struct.atomProperty.macromolecular.residueKey(),
          }),
        data,
      );
      const loci = StructureSelection.toLociWithSourceUnits(sel);
      setStructureOverpaint(
        this.plugin,
        this.plugin.managers.structure.hierarchy.current.structures.flatMap(
          (s) => s.components,
        ),
        Color(color),
        async (structure) => loci,
      );
      this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
    },
    clearHighlight: () => {
      this.plugin.managers.interactivity.lociHighlights.highlightOnly({
        loci: EmptyLoci,
      });
    },
  };
}

(window as any).LouisMolStarWrapper = new BasicWrapper();

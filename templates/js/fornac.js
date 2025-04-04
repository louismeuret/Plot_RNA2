function FornaContainer(h, f) {
  function k() {
    var b = $(h).width(),
      e = $(h).height();
    d.options.svgW = b;
    d.options.svgH = e;
    D.range([0, b]).domain([0, b]);
    E.range([0, e]).domain([0, e]);
    zoomer.x(D).y(E);
    d.brusher.x(D).y(E);
    A.attr("width", b).attr("height", e);
    u.attr("width", b).attr("height", e);
    d.center_view();
  }
  function b(b, d, e) {
    return b.hasOwnProperty(d.num)
      ? ((val = parseFloat(b[d.num])), isNaN(val) ? b[d.num] : e(val))
      : "white";
  }
  function e(b) {
    b = p.selectAll("g.gnode");
    return b.filter(function (b) {
      return b.selected;
    });
  }
  function g(b) {
    var d = b.radius + 16,
      e = b.x - d,
      g = b.x + d,
      f = b.y - d,
      n = b.y + d;
    return function (d, I, l, k, h) {
      if (d.point && d.point !== b) {
        var m = b.x - d.point.x,
          w = b.y - d.point.y,
          q = Math.sqrt(m * m + w * w),
          p = b.radius + d.point.radius;
        q < p &&
          ((q = ((q - p) / q) * 0.1),
          (b.x -= m *= q),
          (b.y -= w *= q),
          (d.point.x += m),
          (d.point.y += w));
      }
      return I > g || k < e || l > n || h < f;
    };
  }
  function n() {
    if (!d.deaf && !C) {
      key_is_down = !0;
      switch (d3.event.keyCode) {
        case 16:
          C = !0;
          break;
        case 17:
          B = !0;
          break;
        case 67:
          d.center_view();
      }
      if (C || B)
        x
          .call(zoomer)
          .on("mousedown.zoom", null)
          .on("touchstart.zoom", null)
          .on("touchmove.zoom", null)
          .on("touchend.zoom", null),
          v.selectAll("g.gnode").on("mousedown.drag", null);
      B &&
        (y.select(".background").style("cursor", "crosshair"),
        y.call(d.brusher));
    }
  }
  function l() {
    B = C = !1;
    y.call(d.brusher)
      .on("mousedown.brush", null)
      .on("touchstart.brush", null)
      .on("touchmove.brush", null)
      .on("touchend.brush", null);
    y.select(".background").style("cursor", "auto");
    x.call(zoomer);
    v.selectAll("g.gnode").call(G);
  }
  var d = this;
  d.options = {
    displayAllLinks: !1,
    labelInterval: 10,
    applyForce: !0,
    initialSize: [200, 200],
    allowPanningAndZooming: !0,
  };
  if (1 < arguments.length)
    for (var w in f) d.options.hasOwnProperty(w) && (d.options[w] = f[w]);
  d.options.svgW = d.options.initialSize[0];
  d.options.svgH = d.options.initialSize[1];
  d3.scale.category20();
  var m = null,
    q = null,
    D = d3.scale
      .linear()
      .domain([0, d.options.svgW])
      .range([0, d.options.svgW]),
    E = d3.scale
      .linear()
      .domain([0, d.options.svgH])
      .range([0, d.options.svgH]),
    F = (d.graph = { nodes: [], links: [] });
  d.linkStrengths = {
    pseudoknot: 0,
    protein_chain: 0,
    chain_chain: 0,
    intermolecule: 10,
    other: 10,
  };
  d.displayParameters = {
    displayBackground: "true",
    displayNumbering: "true",
    displayNodeOutline: "true",
    displayNodeLabel: "true",
    displayLinks: "true",
    displayPseudoknotLinks: "true",
    displayProteinLinks: "true",
  };
  d.colorScheme = "structure";
  d.customColors = {};
  d.animation = d.options.applyForce;
  d.deaf = !1;
  d.rnas = {};
  d.extraLinks = [];
  d.createInitialLayout = function (b, d) {
    var e = {
      sequence: "",
      name: "empty",
      positions: [],
      labelInterval: 10,
      avoidOthers: !0,
    };
    if (2 == arguments.length)
      for (var g in d) e.hasOwnProperty(g) && (e[g] = d[g]);
    rg = new RNAGraph(e.sequence, b, e.name);
    rnaJson = rg.recalculateElements();
    0 === e.positions.length &&
      (e.positions = simple_xy_coordinates(rnaJson.pairtable));
    return (rnaJson = rnaJson
      .elementsToJson()
      .addPositions("nucleotide", e.positions)
      .addLabels(1, e.labelInterval)
      .reinforceStems()
      .reinforceLoops()
      .connectFakeNodes());
  };
  d.addRNA = function (b, e) {
    var g = d.createInitialLayout(b, e);
    1 === arguments.length && (e = {});
    "avoidOthers" in e ? d.addRNAJSON(g, e.avoidOthers) : d.addRNAJSON(g, !0);
    return g;
  };
  d.addRNAJSON = function (b, e) {
    var g, f;
    e &&
      ((g =
        0 < d.graph.nodes.length
          ? d3.max(
              d.graph.nodes.map(function (b) {
                return b.x;
              }),
            )
          : 0),
      (f = d3.min(
        b.nodes.map(function (b) {
          return b.x;
        }),
      )),
      b.nodes.forEach(function (b) {
        b.x += g - f;
        b.px += g - f;
      }));
    b.nodes.forEach(function (d) {
      d.rna = b;
    });
    d.rnas[b.uid] = b;
    d.recalculateGraph();
    d.update();
    d.center_view();
  };
  d.transitionRNA = function (b, e, g) {
    b = d.createInitialLayout(e, g);
    console.log("newRNAJson:", b);
    p.selectAll("g.gnode").each(function (b) {
      console.log("d before", b);
    });
    b = p.selectAll("g.gnode").data(b);
    b.each(function (b) {
      console.log("d after", b);
    });
    b.transition()
      .attr("transform", function (b) {
        console.log("d after", b);
        return "translate(" + [b.x, b.y] + ")";
      })
      .duration(1e3);
  };
  d.recalculateGraph = function (b) {
    d.graph.nodes = [];
    d.graph.links = [];
    for (var e in d.rnas)
      (d.graph.nodes = d.graph.nodes.concat(d.rnas[e].nodes)),
        (d.graph.links = d.graph.links.concat(d.rnas[e].links));
    for (var g = {}, f = 0; f < d.graph.nodes.length; f++)
      g[d.graph.nodes[f].uid] = d.graph.nodes[f];
    d.graph.links.forEach(function (b) {
      b.source = g[b.source.uid];
      b.target = g[b.target.uid];
    });
    for (f = 0; f < d.extraLinks.length; f++) {
      d.extraLinks[f].target.uid in g ||
        console.log("not there:", d.extraLinks[f]);
      d.extraLinks[f].source = g[d.extraLinks[f].source.uid];
      d.extraLinks[f].target = g[d.extraLinks[f].target.uid];
      if ("intermolecule" == d.extraLinks[f].link_type)
        for (
          fake_links = d.graph.links.filter(function (b) {
            return (
              (b.source == d.extraLinks[f].source ||
                b.source == d.extraLinks[f].target ||
                b.target == d.extraLinks[f].source ||
                b.target == d.extraLinks[f].source) &&
              "fake" == b.link_type
            );
          }),
            b = 0;
          b < fake_links.length;
          b++
        )
          (e = d.graph.links.indexOf(fake_links[b])),
            d.graph.links.splice(e, 1);
      F.links.push(d.extraLinks[f]);
    }
  };
  d.addNodes = function (b) {
    b.links.forEach(function (d) {
      "number" == typeof d.source && (d.source = b.nodes[d.source]);
      "number" == typeof d.target && (d.target = b.nodes[d.target]);
    });
    0 < d.graph.nodes.length
      ? ((max_x = d3.max(
          d.graph.nodes.map(function (b) {
            return b.x;
          }),
        )),
        (max_y = d3.max(
          d.graph.nodes.map(function (b) {
            return b.y;
          }),
        )))
      : (max_y = max_x = 0);
    b.nodes.forEach(function (b) {
      b.rna.uid in d.rnas || (d.rnas[b.rna.uid] = b.rna);
      b.x += max_x;
      b.px += max_x;
    });
    r = new RNAGraph("", "");
    r.nodes = b.nodes;
    r.links = b.links;
    d.recalculateGraph();
    d.update();
    d.center_view();
  };
  d.addCustomColors = function (b) {
    d.customColors = b;
  };
  d.clearNodes = function () {
    d.graph.nodes = [];
    d.graph.links = [];
    d.rnas = {};
    d.extraLinks = [];
    d.update();
  };
  d.toJSON = function () {
    return JSON.stringify(
      { rnas: d.rnas, extraLinks: d.extraLinks },
      function (b, d) {
        if ("rna" != b) return d;
      },
      "\t",
    );
  };
  d.fromJSON = function (b) {
    try {
      var e = JSON.parse(b),
        g = e.rnas,
        f = e.extraLinks;
    } catch (n) {
      throw n;
    }
    for (uid in g)
      "rna" == g[uid].type
        ? ((r = new RNAGraph()),
          (r.seq = g[uid].seq),
          (r.dotbracket = g[uid].dotbracket),
          (r.circular = g[uid].circular),
          (r.pairtable = g[uid].pairtable),
          (r.uid = g[uid].uid),
          (r.struct_name = g[uid].struct_name),
          (r.nodes = g[uid].nodes),
          (r.links = g[uid].links),
          (r.rnaLength = g[uid].rnaLength),
          (r.elements = g[uid].elements),
          (r.nucs_to_nodes = g[uid].nucs_to_nodes),
          (r.pseudoknot_pairs = g[uid].pseudoknot_pairs))
        : ((r = new ProteinGraph()),
          (r.size = g[uid].size),
          (r.nodes = g[uid].nodes),
          (r.uid = g[uid].uid)),
        d.addRNAJSON(r, !1);
    f.forEach(function (b) {
      d.extraLinks.push(b);
    });
    d.recalculateGraph();
    d.update();
  };
  d.changeColorScheme = function (e) {
    p.selectAll("[node_type=protein]")
      .classed("protein", !0)
      .attr("r", function (b) {
        return b.radius;
      });
    p.selectAll("g.gnode");
    p.selectAll("g.gnode").selectAll("circle");
    var g = p.selectAll("g.gnode").select("[node_type=nucleotide]");
    d.colorScheme = e;
    "sequence" == e
      ? ((scale = d3.scale
          .ordinal()
          .range(["#dbdb8d", "#98df8a", "#ff9896", "#aec7e8", "#aec7e8"])
          .domain(["A", "C", "G", "U", "T"])),
        g.style("fill", function (b) {
          return scale(b.name);
        }))
      : "structure" == e
        ? ((scale = d3.scale
            .category10()
            .domain("smiethx".split(""))
            .range(
              "lightgreen #ff9896 #dbdb8d lightsalmon lightcyan lightblue transparent".split(
                " ",
              ),
            )),
          g.style("fill", function (b) {
            return scale(b.elem_type);
          }))
        : "positions" == e
          ? g.style("fill", function (b) {
              scale = d3.scale
                .linear()
                .range(["#98df8a", "#dbdb8d", "#ff9896"])
                .interpolate(d3.interpolateLab)
                .domain([1, 1 + (b.rna.rnaLength - 1) / 2, b.rna.rnaLength]);
              return scale(b.num);
            })
          : "custom" == e &&
            ((scale = d3.scale
              .linear()
              .interpolate(d3.interpolateLab)
              .domain(d.customColors.domain)
              .range(d.customColors.range)),
            g.style("fill", function (e) {
              return "undefined" == typeof d.customColors
                ? "white"
                : d.customColors.color_values.hasOwnProperty(e.struct_name) &&
                    d.customColors.color_values[e.struct_name].hasOwnProperty(
                      e.num,
                    )
                  ? ((molecule_colors =
                      d.customColors.color_values[e.struct_name]),
                    b(molecule_colors, e, scale))
                  : d.customColors.color_values.hasOwnProperty("")
                    ? ((molecule_colors = d.customColors.color_values[""]),
                      b(molecule_colors, e, scale))
                    : "white";
            }));
  };
  window.addEventListener("resize", k, !1);
  zoomer = d3.behavior
    .zoom()
    .scaleExtent([0.1, 10])
    .x(D)
    .y(E)
    .on("zoomstart", function () {
      var b = p.selectAll("g.gnode").selectAll(".outline_node");
      b.each(function (b) {
        b.selected = !1;
        b.previouslySelected = !1;
      });
      b.classed("selected", !1);
    })
    .on("zoom", function () {
      v.attr(
        "transform",
        "translate(" + d3.event.translate + ") scale(" + d3.event.scale + ")",
      );
    });
  d3.select(h).select("svg").remove();
  var u = d3
      .select(h)
      .attr("tabindex", 1)
      .on("keydown.brush", n)
      .on("keyup.brush", l)
      .each(function () {
        this.focus();
      })
      .append("svg:svg")
      .attr("width", d.options.svgW)
      .attr("height", d.options.svgH)
      .attr("id", "plotting-area"),
    t = u.append("svg:style");
  $.get(
    location.origin +
      window.location.pathname.replace(/[^\/]+$/g, "") +
      "css/fornac.css",
    function (b) {
      t.text(b.replace(/[\s\n]/g, ""));
    },
  );
  d.options.svg = u;
  var x = u
    .append("svg:g")
    .on("mousemove", function () {
      m &&
        ((mpos = d3.mouse(v.node())),
        H.attr("x1", m.x)
          .attr("y1", m.y)
          .attr("x2", mpos[0])
          .attr("y2", mpos[1]));
    })
    .on("mousedown", function () {})
    .on("mouseup", function () {
      m && H.attr("class", "drag_line_hidden");
      q = m = null;
    });
  d.options.allowPanningAndZooming && x.call(zoomer);
  var A = x
      .append("svg:rect")
      .attr("width", d.options.svgW)
      .attr("height", d.options.svgH)
      .attr("fill", "white")
      .attr("stroke", "transparent")
      .attr("stroke-width", 1)
      .attr("id", "zrect"),
    y = x
      .append("g")
      .datum(function () {
        return { selected: !1, previouslySelected: !1 };
      })
      .attr("class", "brush"),
    v = x.append("svg:g"),
    z = v.append("svg:g"),
    p = v.append("svg:g");
  d.brusher = d3.svg
    .brush()
    .x(D)
    .y(E)
    .on("brushstart", function (b) {
      p.selectAll("g.gnode")
        .selectAll(".outline_node")
        .each(function (b) {
          b.previouslySelected = B && b.selected;
        });
    })
    .on("brush", function () {
      var b = p.selectAll("g.gnode").selectAll(".outline_node"),
        e = d3.event.target.extent();
      b.classed("selected", function (b) {
        return (b.selected =
          d.options.applyForce &&
          b.previouslySelected ^
            (e[0][0] <= b.x &&
              b.x < e[1][0] &&
              e[0][1] <= b.y &&
              b.y < e[1][1]));
      });
    })
    .on("brushend", function () {
      d3.event.target.clear();
      d3.select(this).call(d3.event.target);
    });
  y.call(d.brusher)
    .on("mousedown.brush", null)
    .on("touchstart.brush", null)
    .on("touchmove.brush", null)
    .on("touchend.brush", null);
  y.select(".background").style("cursor", "auto");
  d.center_view = function () {
    0 !== d.graph.nodes.length &&
      ((min_x = d3.min(
        d.graph.nodes.map(function (b) {
          return b.x;
        }),
      )),
      (min_y = d3.min(
        d.graph.nodes.map(function (b) {
          return b.y;
        }),
      )),
      (max_x = d3.max(
        d.graph.nodes.map(function (b) {
          return b.x;
        }),
      )),
      (max_y = d3.max(
        d.graph.nodes.map(function (b) {
          return b.y;
        }),
      )),
      (mol_width = max_x - min_x),
      (mol_height = max_y - min_y),
      (width_ratio = d.options.svgW / (mol_width + 1)),
      (height_ratio = d.options.svgH / (mol_height + 1)),
      (min_ratio = 0.8 * Math.min(width_ratio, height_ratio)),
      (new_mol_width = mol_width * min_ratio),
      (new_mol_height = mol_height * min_ratio),
      (x_trans = -min_x * min_ratio + (d.options.svgW - new_mol_width) / 2),
      (y_trans = -min_y * min_ratio + (d.options.svgH - new_mol_height) / 2),
      v.attr(
        "transform",
        "translate(" + [x_trans, y_trans] + ") scale(" + min_ratio + ")",
      ),
      zoomer.translate([x_trans, y_trans]),
      zoomer.scale(min_ratio));
  };
  d.force = d3.layout
    .force()
    .charge(function (b) {
      return -30;
    })
    .chargeDistance(300)
    .friction(0.35)
    .linkDistance(function (b) {
      return 15 * b.value;
    })
    .linkStrength(function (b) {
      return b.link_type in d.linkStrengths
        ? d.linkStrengths[b.link_type]
        : d.linkStrengths.other;
    })
    .gravity(0)
    .nodes(d.graph.nodes)
    .links(d.graph.links)
    .chargeDistance(110)
    .size([d.options.svgW, d.options.svgH]);
  var H = v
      .append("line")
      .attr("class", "drag_line")
      .attr("x1", 0)
      .attr("y1", 0)
      .attr("x2", 0)
      .attr("y2", 0),
    C = !1,
    B = !1;
  d.resumeForce = function () {
    d.animation && d.force.resume();
  };
  var G = d3.behavior
    .drag()
    .on("dragstart", function (b) {
      d3.event.sourceEvent.stopPropagation();
      b.selected ||
        B ||
        p
          .selectAll("g.gnode")
          .selectAll(".outline_node")
          .classed("selected", function (b) {
            return (b.selected =
              d.options.applyForce && (b.previouslySelected = !1));
          });
      d3.select(this)
        .select(".outline_node")
        .classed("selected", function (e) {
          b.previouslySelected = b.selected;
          return (b.selected = d.options.applyForce && !0);
        });
      e(b).each(function (b) {
        b.fixed |= 2;
      });
    })
    .on("drag", function (b) {
      e(b).each(function (b) {
        b.x += d3.event.dx;
        b.y += d3.event.dy;
        b.px += d3.event.dx;
        b.py += d3.event.dy;
      });
      d.resumeForce();
      d3.event.sourceEvent.preventDefault();
    })
    .on("dragend", function (b) {
      e(b).each(function (b) {
        b.fixed &= -7;
      });
    });
  d3.select(h)
    .on("keydown", n)
    .on("keyup", l)
    .on("contextmenu", function () {
      d3.event.preventDefault();
    });
  link_key = function (b) {
    return b.uid;
  };
  node_key = function (b) {
    return (key = b.uid);
  };
  update_rna_graph = function (b) {
    var e = b.get_positions("nucleotide"),
      g = b.get_positions("label"),
      f = b.get_uids();
    b.recalculateElements()
      .elementsToJson()
      .addPseudoknots()
      .addPositions("nucleotide", e)
      .addUids(f)
      .addLabels(1, d.options.labelInterval)
      .addPositions("label", g)
      .reinforceStems()
      .reinforceLoops()
      .updateLinkUids();
  };
  remove_link = function (b) {
    index = d.graph.links.indexOf(b);
    if (-1 < index) {
      if (b.source.rna == b.target.rna) {
        var e = b.source.rna;
        e.addPseudoknots();
        e.pairtable[b.source.num] = 0;
        e.pairtable[b.target.num] = 0;
        update_rna_graph(e);
      } else
        (extraLinkIndex = d.extraLinks.indexOf(b)),
          d.extraLinks.splice(extraLinkIndex, 1);
      d.recalculateGraph();
    }
    d.update();
  };
  link_click = function (b) {
    C &&
      (b.link_type in
        { backbone: !0, fake: !0, fake_fake: !0, label_link: !0 } ||
        remove_link(b));
  };
  d.add_link = function (b) {
    b.source.rna == b.target.rna
      ? ((r = b.source.rna),
        (r.pairtable[b.source.num] = b.target.num),
        (r.pairtable[b.target.num] = b.source.num),
        update_rna_graph(r))
      : ((b.link_type = "intermolecule"), d.extraLinks.push(b));
    d.recalculateGraph();
    d.update();
  };
  node_mouseclick = function (b) {
    d3.event.defaultPrevented ||
      (B ||
        p
          .selectAll("g.gnode")
          .selectAll(".outline_node")
          .classed("selected", function (b) {
            return (b.selected =
              d.options.applyForce && (b.previouslySelected = !1));
          }),
      d3
        .select(this)
        .select("circle")
        .classed(
          "selected",
          (b.selected = d.options.applyForce && !b.previouslySelected),
        ));
  };
  node_mouseup = function (b) {
    if (m)
      if (((q = b), q == m)) q = m = null;
      else {
        b = {
          source: m,
          target: q,
          link_type: "basepair",
          value: 1,
          uid: generateUUID(),
        };
        for (i = 0; i < d.graph.links.length; i++) {
          if (
            d.graph.links[i].source == m ||
            d.graph.links[i].target == m ||
            d.graph.links[i].source == q ||
            d.graph.links[i].target == q
          )
            if (
              "basepair" == d.graph.links[i].link_type ||
              "pseudoknot" == d.graph.links[i].link_type
            )
              return;
          if (
            ((d.graph.links[i].source == q && d.graph.links[i].target == m) ||
              (d.graph.links[i].source == m && d.graph.links[i].target == q)) &&
            "backbone" == d.graph.links[i].link_type
          )
            return;
        }
        "middle" != q.node_type &&
          "middle" != m.node_type &&
          "label" != q.node_type &&
          "label" != m.node_type &&
          d.add_link(b);
      }
  };
  node_mousedown = function (b) {
    b.selected ||
      B ||
      p
        .selectAll("g.gnode")
        .selectAll(".outline_node")
        .classed("selected", function (b) {
          return (b.selected = b.previouslySelected = !1);
        });
    d3.select(this).classed("selected", function (e) {
      b.previouslySelected = b.selected;
      return (b.selected = d.options.applyForce && !0);
    });
    C &&
      ((m = b),
      H.attr("class", "drag_line")
        .attr("x1", m.x)
        .attr("y1", m.y)
        .attr("x2", m.x)
        .attr("y2", m.y));
  };
  d.startAnimation = function () {
    d.animation = !0;
    v.selectAll("g.gnode").call(G);
    d.force.start();
  };
  d.stopAnimation = function () {
    d.animation = !1;
    v.selectAll("g.gnode").on("mousedown.drag", null);
    d.force.stop();
  };
  d.setFriction = function (b) {
    d.force.friction(b);
    d.resumeForce();
  };
  d.setCharge = function (b) {
    d.force.charge(b);
    d.resumeForce();
  };
  d.setGravity = function (b) {
    d.force.gravity(b);
    d.resumeForce();
  };
  d.setPseudoknotStrength = function (b) {
    d.linkStrengths.pseudoknot = b;
    d.update();
  };
  d.displayBackground = function (b) {
    d.displayParameters.displayBackground = b;
    d.updateStyle();
  };
  d.displayNumbering = function (b) {
    d.displayParameters.displayNumbering = b;
    d.updateStyle();
  };
  d.displayNodeOutline = function (b) {
    d.displayParameters.displayNodeOutline = b;
    d.updateStyle();
  };
  d.displayNodeLabel = function (b) {
    d.displayParameters.displayNodeLabel = b;
    d.updateStyle();
  };
  d.displayLinks = function (b) {
    d.displayParameters.displayLinks = b;
    d.updateStyle();
  };
  d.displayPseudoknotLinks = function (b) {
    d.displayParameters.displayPseudoknotLinks = b;
    d.updateStyle();
  };
  d.displayProteinLinks = function (b) {
    d.displayParameters.displayProteinLinks = b;
    d.updateStyle();
  };
  d.updateStyle = function () {
    A.classed("transparent", !d.displayParameters.displayBackground);
    p.selectAll("[node_type=label]").classed(
      "transparent",
      !d.displayParameters.displayNumbering,
    );
    p.selectAll("[label_type=label]").classed(
      "transparent",
      !d.displayParameters.displayNumbering,
    );
    z.selectAll("[link_type=label_link]").classed(
      "transparent",
      !d.displayParameters.displayNumbering,
    );
    u.selectAll("circle").classed(
      "hidden_outline",
      !d.displayParameters.displayNodeOutline,
    );
    p.selectAll("[label_type=nucleotide]").classed(
      "transparent",
      !d.displayParameters.displayNodeLabel,
    );
    u.selectAll(
      "[link_type=real],[link_type=basepair],[link_type=backbone],[link_type=pseudoknot],[link_type=protein_chain],[link_type=chain_chain]",
    ).classed("transparent", !d.displayParameters.displayLinks);
    u.selectAll("[link_type=pseudoknot]").classed(
      "transparent",
      !d.displayParameters.displayPseudoknotLinks,
    );
    u.selectAll("[link_type=protein_chain]").classed(
      "transparent",
      !d.displayParameters.displayProteinLinks,
    );
    z.selectAll("[link_type=fake]").classed(
      "transparent",
      !d.options.displayAllLinks,
    );
  };
  d.update = function () {
    d.force.nodes(d.graph.nodes).links(d.graph.links);
    d.animation && d.force.start();
    var b = z.selectAll("line.link").data(d.graph.links, link_key);
    link_lines = b.enter().append("svg:line");
    link_lines.append("svg:title").text(link_key);
    link_lines
      .classed("link", !0)
      .attr("x1", function (b) {
        return b.source.x;
      })
      .attr("y1", function (b) {
        return b.source.y;
      })
      .attr("x2", function (b) {
        return b.target.x;
      })
      .attr("y2", function (b) {
        return b.target.y;
      })
      .attr("link_type", function (b) {
        return b.link_type;
      })
      .attr("class", function (b) {
        return d3.select(this).attr("class") + " " + b.link_type;
      })
      .attr("pointer-events", function (b) {
        return "fake" == b.link_type ? "none" : "all";
      });
    b.attr("class", "")
      .classed("link", !0)
      .attr("link_type", function (b) {
        return b.link_type;
      })
      .attr("class", function (b) {
        return d3.select(this).attr("class") + " " + b.link_type;
      });
    b.exit().remove();
    xlink = d.displayFakeLinks
      ? b
      : z.selectAll(
          "[link_type=real],[link_type=pseudoknot],[link_type=protein_chain],[link_type=chain_chain],[link_type=label_link],[link_type=backbone],[link_type=basepair],[link_type=fake],[link_type=intermolecule]",
        );
    domain = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    d3.scale.category10().domain(domain);
    var e = p.selectAll("g.gnode").data(d.graph.nodes, node_key);
    gnodes_enter = e
      .enter()
      .append("g")
      .classed("noselect", !0)
      .classed("gnode", !0)
      .attr("struct_name", function (b) {
        return b.struct_name;
      })
      .attr("transform", function (b) {
        return "undefined" != typeof b.x && "undefined" != typeof b.y
          ? "translate(" + [b.x, b.y] + ")"
          : "";
      })
      .each(function (b) {
        b.selected = b.previouslySelected = !1;
      });
    gnodes_enter
      .call(G)
      .on("mousedown", node_mousedown)
      .on("mousedrag", function (b) {})
      .on("mouseup", node_mouseup)
      .on("click", node_mouseclick)
      .transition()
      .duration(750)
      .ease("elastic")
      .attr("r", 6.5);
    node_tooltip = function (b) {
      node_tooltips = {};
      node_tooltips.nucleotide = b.num;
      node_tooltips.label = "";
      node_tooltips.pseudo = "";
      node_tooltips.middle = "";
      node_tooltips.protein = b.struct_name;
      return node_tooltips[b.node_type];
    };
    xlink.on("click", link_click);
    e.select("circle");
    gnodes_enter
      .filter(function (b) {
        return (
          "nucleotide" == b.node_type ||
          "label" == b.node_type ||
          "protein" == b.node_type
        );
      })
      .append("svg:circle")
      .attr("class", "outline_node")
      .attr("r", function (b) {
        return b.radius + 1;
      });
    b = gnodes_enter
      .append("svg:circle")
      .attr("class", "node")
      .classed("label", function (b) {
        return "label" == b.node_type;
      })
      .attr("r", function (b) {
        return "middle" == b.node_type ? 0 : b.radius;
      })
      .attr("node_type", function (b) {
        return b.node_type;
      });
    gnodes_enter
      .append("text")
      .text(function (b) {
        return b.name;
      })
      .attr("text-anchor", "middle")
      .attr("font-size", 8)
      .attr("font-weight", "bold")
      .attr("y", 2.5)
      .attr("class", "node-label")
      .attr("label_type", function (b) {
        return b.node_type;
      })
      .append("svg:title")
      .text(function (b) {
        return "nucleotide" == b.node_type ? b.struct_name + ":" + b.num : "";
      });
    b.append("svg:title").text(function (b) {
      return "nucleotide" == b.node_type ? b.struct_name + ":" + b.num : "";
    });
    e.exit().remove();
    real_nodes = d.graph.nodes.filter(function (b) {
      return "nucleotide" == b.node_type || "label" == b.node_type;
    });
    d.force.on("tick", function () {
      for (
        var b = d3.geom.quadtree(real_nodes), d = 0, f = real_nodes.length;
        ++d < f;

      )
        b.visit(g(real_nodes[d]));
      xlink
        .attr("x1", function (b) {
          return b.source.x;
        })
        .attr("y1", function (b) {
          return b.source.y;
        })
        .attr("x2", function (b) {
          return b.target.x;
        })
        .attr("y2", function (b) {
          return b.target.y;
        });
      e.attr("transform", function (b) {
        return "translate(" + [b.x, b.y] + ")";
      });
    });
    d.changeColorScheme(d.colorScheme);
    d.animation && d.force.start();
    d.updateStyle();
  };
  k();
}
var number_sort = function (h, f) {
  return h - f;
};
function generateUUID() {
  var h = new Date().getTime();
  return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g, function (f) {
    var k = (h + 16 * Math.random()) % 16 | 0;
    h = Math.floor(h / 16);
    return ("x" == f ? k : (k & 3) | 8).toString(16);
  });
}
function isNormalInteger(h) {
  return /^\+?(0|[1-9]\d*)$/.test(h);
}
"undefined" === typeof String.prototype.trim &&
  (String.prototype.trim = function () {
    return String(this).replace(/^\s+|\s+$/g, "");
  });
function ColorScheme(h) {
  var f = this;
  f.colorsText = h;
  f.parseRange = function (f) {
    for (var b = f.split(","), e = [], g = 0; g < b.length; g++) {
      var n = b[g].split("-");
      if (1 == n.length) e.push(parseInt(n[0]));
      else if (2 == n.length)
        for (var l = parseInt(n[0]), n = parseInt(n[1]); l <= n; l++) e.push(l);
      else console.log("Malformed range (too many dashes):", f);
    }
    return e;
  };
  f.parseColorText = function (h) {
    h = h.split("\n");
    for (
      var b = "",
        e = 1,
        g = { color_values: { "": {} }, range: ["white", "steelblue"] },
        n = [],
        l = 0;
      l < h.length;
      l++
    )
      if (">" == h[l][0])
        (b = h[l].trim().slice(1)), (e = 1), (g.color_values[b] = {});
      else {
        words = h[l].trim().split(/[\s]+/);
        for (var d = 0; d < words.length; d++)
          if (isNaN(words[d]))
            if (0 === words[d].search("range"))
              (parts = words[d].split("=")),
                (parts_right = parts[1].split(":")),
                (g.range = [parts_right[0], parts_right[1]]);
            else if (0 == words[d].search("domain"))
              (parts = words[d].split("=")),
                (parts_right = parts[1].split(":")),
                (g.domain = [parts_right[0], parts_right[1]]);
            else {
              parts = words[d].split(":");
              nums = f.parseRange(parts[0]);
              color = parts[1];
              for (var w = 0; w < nums.length; w++)
                isNaN(color)
                  ? (g.color_values[b][nums[w]] = color)
                  : ((g.color_values[b][nums[w]] = +color),
                    n.push(Number(color)));
            }
          else
            (g.color_values[b][e] = Number(words[d])),
              (e += 1),
              n.push(Number(words[d]));
      }
    "domain" in g ||
      (g.domain = [Math.min.apply(null, n), Math.max.apply(null, n)]);
    f.colors_json = g;
    return f;
  };
  f.normalizeColors = function () {
    var h, b;
    for (b in f.colors_json) {
      var e = Number.MAX_VALUE,
        g = Number.MIN_VALUE,
        n;
      for (n in f.colors_json.color_values[b])
        (h = f.colors_json.color_values[b][n]),
          "number" == typeof h && (h < e && (e = h), h > g && (g = h));
      for (n in f.colors_json.color_values[b])
        (h = f.colors_json.color_values[b][n]),
          "number" == typeof h &&
            (f.colors_json.color_values[b][n] = (h - e) / (g - e));
    }
    return f;
  };
  f.parseColorText(f.colorsText);
}
function ProteinGraph(h, f, k) {
  var b = this;
  b.type = "protein";
  b.size = f;
  b.nodes = [
    {
      name: "P",
      num: 1,
      radius: 3 * Math.sqrt(f),
      rna: b,
      node_type: "protein",
      struct_name: h,
      elem_type: "p",
      size: f,
      uid: generateUUID(),
    },
  ];
  b.links = [];
  b.uid = generateUUID();
  b.addUids = function (e) {
    for (var g = 0; g < e.length; g++) b.nodes[g].uid = e[g];
    return b;
  };
  b.get_uids = function () {
    uids = [];
    for (var e = 0; e < b.dotbracket.length; e++) uids.push(b.nodes[e].uid);
    return uids;
  };
}
function RNAGraph(h, f, k) {
  var b = this;
  b.type = "rna";
  b.circularizeExternal = !1;
  0 == arguments.length
    ? ((b.seq = ""), (b.dotbracket = ""), (b.struct_name = ""))
    : ((b.seq = h), (b.dotbracket = f), (b.struct_name = k));
  b.circular = !1;
  0 < b.dotbracket.length &&
    "*" == b.dotbracket[b.dotbracket.length - 1] &&
    ((b.dotbracket = b.dotbracket.slice(0, b.dotbracket.length - 1)),
    (b.circular = !0));
  b.uid = generateUUID();
  b.rnaLength = b.dotbracket.length;
  b.elements = [];
  b.pseudoknotPairs = [];
  b.nucs_to_nodes = {};
  b.addUids = function (e) {
    for (var g = 0; g < e.length; g++) b.nodes[g].uid = e[g];
    return b;
  };
  b.computePairtable = function () {
    b.pairtable = rnaUtilities.dotbracketToPairtable(b.dotbracket);
  };
  b.computePairtable();
  b.addPositions = function (e, g) {
    label_nodes = b.nodes.filter(function (b) {
      return b.node_type == e;
    });
    for (var f = 0; f < label_nodes.length; f++)
      (label_nodes[f].x = g[f][0]),
        (label_nodes[f].px = g[f][0]),
        (label_nodes[f].y = g[f][1]),
        (label_nodes[f].py = g[f][1]);
    return b;
  };
  b.get_positions = function (e) {
    positions = [];
    nucleotide_nodes = b.nodes.filter(function (b) {
      return b.node_type == e;
    });
    for (var g = 0; g < nucleotide_nodes.length; g++)
      positions.push([nucleotide_nodes[g].x, nucleotide_nodes[g].y]);
    return positions;
  };
  b.get_uids = function () {
    uids = [];
    for (var e = 0; e < b.dotbracket.length; e++) uids.push(b.nodes[e].uid);
    return uids;
  };
  b.reinforceStems = function () {
    pt = b.pairtable;
    relevant_elements = elements.filter(function (b) {
      return "s" == b[0] && 4 <= b[2].length;
    });
    for (var e = 0; e < relevant_elements.length; e++) {
      all_nucs = relevant_elements[e][2];
      nucs = all_nucs.slice(0, all_nucs.length / 2);
      for (var g = 0; g < nucs.length - 1; g++)
        b.addFakeNode([nucs[g], nucs[g + 1], pt[nucs[g + 1]], pt[nucs[g]]]);
    }
    return b;
  };
  b.reinforceLoops = function () {
    var e = function (d) {
      return 0 !== d && d <= b.dotbracket.length;
    };
    for (i = 0; i < b.elements.length; i++)
      if (
        "s" != b.elements[i][0] &&
        (b.circularizeExternal || "e" != b.elements[i][0])
      ) {
        var g = b.elements[i][2].filter(e);
        if ("e" == b.elements[i][0]) {
          var f = {
              name: "",
              num: -1,
              radius: 0,
              rna: b,
              node_type: "middle",
              elem_type: "f",
              nucs: [],
              x: b.nodes[b.rnaLength - 1].x,
              y: b.nodes[b.rnaLength - 1].y,
              px: b.nodes[b.rnaLength - 1].px,
              py: b.nodes[b.rnaLength - 1].py,
              uid: generateUUID(),
            },
            h = {
              name: "",
              num: -1,
              radius: 0,
              rna: b,
              node_type: "middle",
              elem_type: "f",
              nucs: [],
              x: b.nodes[0].x,
              y: b.nodes[0].y,
              px: b.nodes[0].px,
              py: b.nodes[0].py,
              uid: generateUUID(),
            };
          g.push(b.nodes.length + 1);
          g.push(b.nodes.length + 2);
          b.nodes.push(f);
          b.nodes.push(h);
        }
        b.addFakeNode(g);
      }
    return b;
  };
  b.updateLinkUids = function () {
    for (var e = 0; e < b.links.length; e++)
      b.links[e].uid = b.links[e].source.uid + b.links[e].target.uid;
    return b;
  };
  b.addFakeNode = function (e) {
    for (
      var g = 6.283 / (2 * e.length), g = 18 / (2 * Math.tan(g)), f = "", h = 0;
      h < e.length;
      h++
    )
      f += b.nodes[e[h] - 1].uid;
    new_node = {
      name: "",
      num: -1,
      radius: g,
      rna: b,
      node_type: "middle",
      elem_type: "f",
      nucs: e,
      uid: f,
    };
    b.nodes.push(new_node);
    coords_counted = new_y = new_x = 0;
    g = (3.14159 * (e.length - 2)) / (2 * e.length);
    g = 0.5 / Math.cos(g);
    for (j = 0; j < e.length; j++)
      0 === e[j] ||
        e[j] > b.dotbracket.length ||
        (b.links.push({
          source: b.nodes[e[j] - 1],
          target: b.nodes[b.nodes.length - 1],
          link_type: "fake",
          value: g,
          uid: generateUUID(),
        }),
        4 < e.length &&
          b.links.push({
            source: b.nodes[e[j] - 1],
            target: b.nodes[e[(j + Math.floor(e.length / 2)) % e.length] - 1],
            link_type: "fake",
            value: 2 * g,
            uid: generateUUID(),
          }),
        (ia = (3.14159 * (e.length - 2)) / e.length),
        (c = 2 * Math.cos(1.570795 - ia / 2)),
        b.links.push({
          source: b.nodes[e[j] - 1],
          target: b.nodes[e[(j + 2) % e.length] - 1],
          link_type: "fake",
          value: c,
        }),
        (from_node = b.nodes[e[j] - 1]),
        "x" in from_node &&
          ((new_x += from_node.x),
          (new_y += from_node.y),
          (coords_counted += 1)));
    0 < coords_counted &&
      ((new_node.x = new_x / coords_counted),
      (new_node.y = new_y / coords_counted),
      (new_node.px = new_node.x),
      (new_node.py = new_node.y));
    return b;
  };
  b.connectFakeNodes = function () {
    for (
      var e = {},
        g = b.nodes.filter(function (b) {
          return "middle" == b.node_type;
        }),
        f = new Set(),
        h = 1;
      h <= b.nodes.length;
      h++
    )
      e[h] = [];
    for (h = 0; h < g.length; h++)
      for (var d = g[h], k = 0; k < d.nucs.length; k++) {
        for (var m = d.nucs[k], q = 0; q < e[m].length; q++)
          f.has(JSON.stringify([e[m][q].uid, d.uid].sort())) ||
            (b.links.push({
              source: e[m][q],
              target: d,
              value: (e[m][q].radius + d.radius) / 18,
              link_type: "fake_fake",
            }),
            f.add(JSON.stringify([e[m][q].uid, d.uid].sort())));
        e[m].push(d);
      }
    return b;
  };
  b.elementsToJson = function () {
    pt = b.pairtable;
    elements = b.elements;
    b.nodes = [];
    b.links = [];
    elem_types = {};
    b.elements.sort();
    for (var e = 0; e < b.elements.length; e++)
      for (nucs = b.elements[e][2], j = 0; j < nucs.length; j++)
        elem_types[nucs[j]] = b.elements[e][0];
    for (e = 1; e <= pt[0]; e++)
      b.nodes.push({
        name: b.seq[e - 1],
        num: e,
        radius: 6,
        rna: b,
        node_type: "nucleotide",
        struct_name: b.struct_name,
        elem_type: elem_types[e],
        uid: generateUUID(),
      });
    for (e = 1; e <= pt[0]; e++)
      0 !== pt[e] &&
        b.links.push({
          source: b.nodes[e - 1],
          target: b.nodes[pt[e] - 1],
          link_type: "basepair",
          value: 1,
          uid: generateUUID(),
        }),
        1 < e &&
          b.links.push({
            source: b.nodes[e - 2],
            target: b.nodes[e - 1],
            link_type: "backbone",
            value: 1,
            uid: generateUUID(),
          });
    for (e = 0; e < b.pseudoknotPairs.length; e++)
      b.links.push({
        source: b.nodes[b.pseudoknotPairs[e][0] - 1],
        target: b.nodes[b.pseudoknotPairs[e][1] - 1],
        link_type: "pseudoknot",
        value: 1,
        uid: generateUUID(),
      });
    b.circular &&
      b.links.push({
        source: b.nodes[0],
        target: b.nodes[b.rnaLength - 1],
        link_type: "backbone",
        value: 1,
        uid: generateUUID(),
      });
    return b;
  };
  b.pt_to_elements = function (e, g, f, h) {
    var d = [],
      k = [f - 1],
      m = [h + 1];
    if (f > h) return [];
    for (; 0 === e[f]; f++) k.push(f);
    for (; 0 === e[h]; h--) m.push(h);
    if (f > h)
      return (
        k.push(f),
        0 === g
          ? [["e", g, k.sort(number_sort)]]
          : [["h", g, k.sort(number_sort)]]
      );
    if (e[f] != h) {
      for (k.push(f); f <= h; ) {
        d = d.concat(b.pt_to_elements(e, g, f, e[f]));
        k.push(e[f]);
        for (f = e[f] + 1; 0 === e[f] && f <= h; f++) k.push(f);
        k.push(f);
      }
      k.pop();
      k = k.concat(m);
      0 < k.length &&
        (0 === g
          ? d.push(["e", g, k.sort(number_sort)])
          : d.push(["m", g, k.sort(number_sort)]));
      return d;
    }
    e[f] === h &&
      (k.push(f),
      m.push(h),
      (combined = k.concat(m)),
      4 < combined.length &&
        (0 === g
          ? d.push(["e", g, k.concat(m).sort(number_sort)])
          : d.push(["i", g, k.concat(m).sort(number_sort)])));
    for (m = []; e[f] === h && f < h; )
      m.push(f), m.push(h), (f += 1), --h, (g += 1);
    d.push(["s", g, m.sort(number_sort)]);
    return d.concat(b.pt_to_elements(e, g, f, h));
  };
  b.addLabels = function (e, g) {
    0 === arguments.length && ((e = 1), (g = 10));
    1 === arguments.length && (g = 10);
    if (0 === g) return b;
    0 >= g && console.log("The label interval entered in invalid:", g);
    for (i = 1; i <= pt[0]; i++)
      if (0 === i % g) {
        var f, h;
        thisNode = b.nodes[i - 1];
        1 == b.rnaLength
          ? ((nextVec = [thisNode.x - 15, thisNode.y]),
            (prevVec = [thisNode.x - 15, thisNode.y]))
          : ((prevNode = 1 == i ? b.nodes[b.rnaLength - 1] : b.nodes[i - 2]),
            (nextNode = i == b.rnaLength ? b.nodes[0] : b.nodes[i]),
            0 !== b.pairtable[nextNode.num] &&
              0 !== b.pairtable[prevNode.num] &&
              0 !== b.pairtable[thisNode.num] &&
              (prevNode = nextNode = b.nodes[b.pairtable[thisNode.num] - 1]),
            0 === b.pairtable[thisNode.num] ||
            (0 !== b.pairtable[nextNode.num] && 0 !== b.pairtable[prevNode.num])
              ? ((nextVec = [nextNode.x - thisNode.x, nextNode.y - thisNode.y]),
                (prevVec = [prevNode.x - thisNode.x, prevNode.y - thisNode.y]))
              : ((nextVec = [thisNode.x - nextNode.x, thisNode.y - nextNode.y]),
                (prevVec = [
                  thisNode.x - prevNode.x,
                  thisNode.y - prevNode.y,
                ])));
        combinedVec = [nextVec[0] + prevVec[0], nextVec[1] + prevVec[1]];
        vecLength = Math.sqrt(
          combinedVec[0] * combinedVec[0] + combinedVec[1] * combinedVec[1],
        );
        normedVec = [combinedVec[0] / vecLength, combinedVec[1] / vecLength];
        offsetVec = [-15 * normedVec[0], -15 * normedVec[1]];
        f = b.nodes[i - 1].x + offsetVec[0];
        h = b.nodes[i - 1].y + offsetVec[1];
        new_node = {
          name: i + e - 1,
          num: -1,
          radius: 6,
          rna: b,
          node_type: "label",
          struct_name: b.struct_name,
          elem_type: "l",
          x: f,
          y: h,
          px: f,
          py: h,
          uid: generateUUID(),
        };
        new_link = {
          source: b.nodes[i - 1],
          target: new_node,
          value: 1,
          link_type: "label_link",
          uid: generateUUID(),
        };
        b.nodes.push(new_node);
        b.links.push(new_link);
      }
    return b;
  };
  b.recalculateElements = function () {
    b.removePseudoknots();
    b.elements = b.pt_to_elements(b.pairtable, 0, 1, b.dotbracket.length);
    if (
      b.circular &&
      ((external_loop = b.elements.filter(function (b) {
        if ("e" == b[0]) return !0;
      })),
      0 < external_loop.length)
    ) {
      eloop = external_loop[0];
      nucs = eloop[2].sort(number_sort);
      prev = nucs[0];
      hloop = !0;
      num_greater = 0;
      for (var e = 1; e < nucs.length; e++)
        1 < nucs[e] - prev && (num_greater += 1), (prev = nucs[e]);
      eloop[0] = 1 == num_greater ? "h" : 2 == num_greater ? "i" : "m";
    }
    return b;
  };
  b.removePseudoknots = function () {
    b.pseudoknotPairs =
      1 < b.pairtable.length
        ? rnaUtilities.removePseudoknotsFromPairtable(b.pairtable)
        : [];
    return b;
  };
  b.addPseudoknots = function () {
    var e = b.pairtable,
      f = b.pseudoknotPairs;
    for (i = 0; i < f.length; i++)
      (e[f[i][0]] = f[i][1]), (e[f[i][1]] = f[i][0]);
    b.pseudoknotPairs = [];
    return b;
  };
  0 < b.rnaLength && b.recalculateElements();
}
molecules_to_json = function (h) {
  for (var f = {}, k = [], b = [], e = 0; e < h.molecules.length; e++) {
    var g = h.molecules[e];
    "rna" == g.type
      ? ((rg = new RNAGraph(g.seq, g.ss, g.header)),
        (rg.circularizeExternal = !0),
        rg
          .elementsToJson()
          .addPositions("nucleotide", g.positions)
          .addLabels()
          .reinforceStems()
          .reinforceLoops())
      : "protein" == g.type && (rg = new ProteinGraph(g.header, g.size));
    rg.addUids(g.uids);
    for (g = 0; g < rg.nodes.length; g++) f[rg.nodes[g].uid] = rg.nodes[g];
    k.push(rg);
  }
  for (e = 0; e < h.extra_links.length; e++)
    (link = h.extra_links[e]),
      (link.source = f[link.source]),
      (link.target = f[link.target]),
      (link.uid = generateUUID()),
      b.push(link);
  return { graphs: k, extraLinks: b };
};
number_sort = function (h, f) {
  return h - f;
};
function RNAUtilities() {
  var h = this;
  h.bracket_left = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
  h.bracket_right = ")]}>abcdefghijklmnopqrstuvwxyz".split("");
  h.inverse_brackets = function (f) {
    res = {};
    for (i = 0; i < f.length; i++) res[f[i]] = i;
    return res;
  };
  h.maximumMatching = function (f) {
    var h = f[0];
    mm = Array(h + 1);
    for (var b = 0; b <= h; b++) {
      mm[b] = Array(h + 1);
      for (var e = b; e <= h; e++) mm[b][e] = 0;
    }
    for (var g = 0, b = h - 0 - 1; 0 < b; b--)
      for (e = b + 0 + 1; e <= h; e++) {
        for (var g = mm[b][e - 1], n = e - 0 - 1; n >= b; n--)
          f[n] === e &&
            (g = Math.max(
              g,
              (n > b ? mm[b][n - 1] : 0) +
                1 +
                (0 < e - n - 1 ? mm[n + 1][e - 1] : 0),
            ));
        mm[b][e] = g;
      }
    return mm;
  };
  h.backtrackMaximumMatching = function (f, k) {
    var b = Array.apply(null, Array(f.length)).map(function () {
      return 0;
    });
    h.mm_bt(f, b, k, 1, f.length - 1);
    return b;
  };
  h.mm_bt = function (f, k, b, e, g) {
    var n = f[e][g];
    if (!(0 > g - e - 1))
      if (f[e][g - 1] == n) h.mm_bt(f, k, b, e, g - 1);
      else {
        for (var l = g - 0 - 1; l >= e; l--)
          if (
            b[g] === l &&
            (l > e ? f[e][l - 1] : 0) +
              (0 < g - l - 1 ? f[l + 1][g - 1] : 0) +
              1 ==
              n
          ) {
            k[l] = g;
            k[g] = l;
            e < l && h.mm_bt(f, k, b, e, l - 1);
            h.mm_bt(f, k, b, l + 1, g - 1);
            return;
          }
        console.log("FAILED!!!" + e + "," + g + ": backtracking failed!");
      }
  };
  h.dotbracketToPairtable = function (f) {
    pt = Array.apply(null, Array(f.length + 1)).map(
      Number.prototype.valueOf,
      0,
    );
    pt[0] = f.length;
    stack = {};
    for (i = 0; i < h.bracket_left.length; i++) stack[i] = [];
    inverse_bracket_left = h.inverse_brackets(h.bracket_left);
    inverse_bracket_right = h.inverse_brackets(h.bracket_right);
    for (i = 0; i < f.length; i++)
      if (((a = f[i]), (ni = i + 1), "." == a)) pt[ni] = 0;
      else if (a in inverse_bracket_left)
        stack[inverse_bracket_left[a]].push(ni);
      else if (a in inverse_bracket_right)
        (j = stack[inverse_bracket_right[a]].pop()), (pt[ni] = j), (pt[j] = ni);
      else throw "Unknown symbol in dotbracket string";
    for (key in stack)
      if (0 < stack[key].length)
        throw "Unmatched base at position " + stack[key][0];
    return pt;
  };
  h.insert_into_stack = function (f, h, b) {
    for (h = 0; 0 < f[h].length && f[h][f[h].length - 1] < b; ) h += 1;
    f[h].push(b);
    return h;
  };
  h.delete_from_stack = function (f, h) {
    for (var b = 0; 0 === f[b].length || f[b][f[b].length - 1] != h; ) b += 1;
    f[b].pop();
    return b;
  };
  h.pairtableToDotbracket = function (f) {
    stack = {};
    for (i = 0; i < f[0]; i++) stack[i] = [];
    seen = {};
    res = "";
    for (i = 1; i < f[0] + 1; i++) {
      if (0 !== f[i] && f[i] in seen)
        throw "Invalid pairtable contains duplicate entries";
      seen[f[i]] = !0;
      res =
        0 === f[i]
          ? res + "."
          : f[i] > i
            ? res + h.bracket_left[h.insert_into_stack(stack, i, f[i])]
            : res + h.bracket_right[h.delete_from_stack(stack, i)];
    }
    return res;
  };
  h.find_unmatched = function (f, k, b) {
    for (var e = [], g = [], n = b, l = k; l <= b; l++)
      0 !== f[l] && (f[l] < k || f[l] > b) && g.push([l, f[l]]);
    for (l = k; l <= n; l++) {
      for (; 0 === f[l] && l <= n; ) l++;
      for (b = f[l]; f[l] === b; ) l++, b--;
      e = e.concat(h.find_unmatched(f, l, b));
    }
    0 < g.length && e.push(g);
    return e;
  };
  h.removePseudoknotsFromPairtable = function (f) {
    for (
      var k = h.maximumMatching(f),
        k = h.backtrackMaximumMatching(k, f),
        b = [],
        e = 1;
      e < f.length;
      e++
    )
      f[e] < e ||
        k[e] == f[e] ||
        (b.push([e, f[e]]), (f[f[e]] = 0), (f[e] = 0));
    return b;
  };
}
rnaUtilities = new RNAUtilities();
simple_xy_coordinates = function (h) {
  var f = [],
    k = [];
  console.log("pair_table", h);
  var b, e;
  b = h[0];
  var g = Array.apply(null, Array(b + 5)).map(Number.prototype.valueOf, 0),
    n = Array.apply(null, Array(16 + Math.floor(b / 5))).map(
      Number.prototype.valueOf,
      0,
    ),
    l = Array.apply(null, Array(16 + Math.floor(b / 5))).map(
      Number.prototype.valueOf,
      0,
    );
  lp = stk = 0;
  var d = Math.PI / 2;
  loop = function (b, e, f) {
    var h = 2,
      k = 0,
      F = 0,
      u,
      t,
      x,
      A,
      y,
      v,
      z;
    console.log("i:", b, "j:", e);
    var p = Array.apply(null, Array(1 + 2 * Math.floor((e - b) / 5))).map(
      Number.prototype.valueOf,
      0,
    );
    u = b - 1;
    for (e++; b != e; )
      if ((t = f[b]) && 0 != b) {
        h += 2;
        x = b;
        A = t;
        p[++k] = x;
        p[++k] = A;
        b = t + 1;
        y = x;
        v = A;
        z = 0;
        do x++, A--, z++;
        while (f[x] == A);
        t = z - 2;
        if (
          2 <= z &&
          ((g[y + 1 + t] += d),
          (g[v - 1 - t] += d),
          (g[y] += d),
          (g[v] += d),
          2 < z)
        )
          for (; 1 <= t; t--) (g[y + t] = Math.PI), (g[v - t] = Math.PI);
        l[++stk] = z;
        loop(x, A, f);
      } else b++, h++, F++;
    b = (Math.PI * (h - 2)) / h;
    p[++k] = e;
    e = 0 > u ? 0 : u;
    for (u = 1; u <= k; u++) {
      f = p[u] - e;
      for (t = 0; t <= f; t++) g[e + t] += b;
      if (u > k) break;
      e = p[++u];
    }
    n[++lp] = F;
  };
  loop(0, b + 1, h);
  n[lp] -= 2;
  e = 0;
  f[0] = 100;
  k[0] = 100;
  poss = [];
  poss.push([f[0], k[0]]);
  for (h = 1; h < b; h++)
    (f[h] = f[h - 1] + 15 * Math.cos(e)),
      (k[h] = k[h - 1] + 15 * Math.sin(e)),
      poss.push([f[h], k[h]]),
      (e += Math.PI - g[h + 1]);
  return poss;
};

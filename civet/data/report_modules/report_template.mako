<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="icon" href="https://raw.githubusercontent.com/cov-ert/civet/master/docs/virus.svg">

    <title>${config["report_title"]}</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha384-HSMxcRTRxnN+Bdg0JdbxYKrThecOKuH5zCYotlSAcp1+c8xmyTe9GYg1l9a69psu" crossorigin="anonymous">
    <script src="https://code.jquery.com/jquery-1.12.4.min.js" integrity="sha384-nvAa0+6Qg9clwYCGGPpDQLVpLNn0fRaROjHqs13t4Ggj3Ez50XnGQqc/r8MhnRDZ" crossorigin="anonymous"></script>
    <link href="https://cdn.datatables.net/1.10.25/css/jquery.dataTables.min.css" rel="stylesheet" type="text/css" />
    <script src="https://cdn.datatables.net/1.10.25/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/dataTables.buttons.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.html5.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.print.min.js"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha384-aJ21OjlMXNL5UyIl/XNwTMqvzeRMZH2w8c5cRVpzpU8Y5bApTppSuUkhZXN0VxHd" crossorigin="anonymous"></script>
    <script src="https://cdn.jsdelivr.net/gh/rambaut/figtree.js@9880/dist/figtree.umd.js"></script>
    <script src="https://d3js.org/d3.v6.min.js"></script>
    <script src="https://sharonchoong.github.io/svg-exportJS/svg-export.min.js"></script>
    <script src="https://unpkg.com/canvg/lib/umd.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/pdfkit/js/pdfkit.min.js"></script>
    <script src="https://github.com/devongovett/blob-stream/releases/download/v0.1.3/blob-stream.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/svg-to-pdfkit/source.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega@5.16.0"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@4.15.0"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@6.11.1"></script>

    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <% colorCodes = config["colour_map"] %>
    <% themeColor = config["colour_theme"] %>
    <style>
      body {
        padding-top: 50px;
        font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif;
      }
      table text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
      }
      header {
          display: block;
          text-align: right;
      
      }
      svg {width: 90%; height: auto;}
      .row {
          display: flex;
        }
        .column {
          padding: 10px;
          flex: 50%;
        }
      .accordion {
          background-color: #eee;
          color: #444;
          cursor: pointer;
          padding: 13px;
          width: 100%;
          border: none;
          text-align: left;
          outline: none;
          transition: 0.4s;
        }

        .active, .accordion:hover {
          background-color: ${themeColor};
          color: white;
        }

        .accordion:after {
          content: '\002B';
          color: white;
          font-weight: bold;
          float: right;
          margin-left: 5px;
        }

        .active:after {
          content: "\2212";
        }

        .panel {
          padding: 0 13px;
          background-color: white;
          max-height: 0;
          overflow: hidden;
          transition: max-height 0.2s ease-out;
        }
      .center {
          display: block;
          margin-left: auto;
          margin-right: auto;
          width: 50%;
          }
      .node-background{
          fill:dimgrey;
          stroke:dimgrey;
      }
      .node circle{
        stroke-width:0;
        cursor:pointer;
        /* fill:#7178bc; */
        stroke:dimgrey;
        }
      .node circle.selected{
        stroke-width:0;
        cursor:pointer;
        fill:${themeColor};
        stroke:dimgrey;
        }
      .node-background.query_boolean-True{
          stroke:${themeColor};
      }
      .node.query_boolean-True circle{
        stroke:${themeColor};
      }
      .node.query_boolean-True circle.selected{
        stroke:${themeColor};
      }
      .node-background.query_boolean-True circle.selected{
          stroke:${themeColor};
      }
      .node.query_boolean-True.hovered circle{
          stroke:${themeColor};
      }
      .node rect{
        stroke-width:2;
        fill:${themeColor};
        stroke:dimgrey;
      }
      .svg-tooltip {
          background: rgba(69,77,93,.9);
          border-radius: .1rem;
          color: #fff;
          display: block;
          font-size: 14px;
          max-width: 320px;
          padding: .2rem .4rem;
          position: absolute;
          text-overflow: ellipsis;
          white-space: pre;
          z-index: 300;
          visibility: hidden;
    }
    .tooltip-header {
      font-size: 1.3em;
    }
    .tooltip-key {
      font-weight: bold;
    }
    .branch path{
      stroke-width:2;
      stroke: dimgrey;
      stroke-linejoin:round;
      cursor: pointer;
      }
      .branch.hovered path{
        stroke-width:4;
        stroke: dimgrey;
      }
        .node.hovered circle{
        stroke-width:5;
        stroke: dimgrey
        }
        .node text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
          font-weight: 300;
          font-size: 0.9em;
        }
      /* .starter-template {
        padding: 40px 15px;
        text-align: left;
      } */
      .dataTables_wrapper.no-footer .dataTables_scrollBody {
        border-top: 1px solid  rgb(148, 148, 148);
        border-bottom: none;
      }
      .svg-icon {
      display: inline-flex;
      align-self: center;
      }
      h3{
          font-size: 1em;
      }
      #toTopBtn {
      position: fixed;
      bottom: 26px;
      right: 39px;
      z-index: 98;
      padding: 21px;
      background-color: ${themeColor}
      }
      .js .cd-top--fade-out {
          opacity: .5
      }
      .js .cd-top--is-visible {
          visibility: visible;
          opacity: 1
      }
      .js .cd-top {
          visibility: hidden;
          opacity: 0;
          transition: opacity .3s, visibility .3s, background-color .3s
      }
      .cd-top {
          position: fixed;
          bottom: 20px;
          bottom: var(--cd-back-to-top-margin);
          right: 20px;
          right: var(--cd-back-to-top-margin);
          display: inline-block;
          height: 40px;
          height: var(--cd-back-to-top-size);
          width: 40px;
          width: var(--cd-back-to-top-size);
          box-shadow: 0 0 10px rgba(0, 0, 0, .05) !important;
          background: url(https://res.cloudinary.com/dxfq3iotg/image/upload/v1571057658/cd-top-arrow.svg) no-repeat center 50%;
          background-color: ${themeColor};
          background-color: hsla(var(--cd-color-3-h), var(--cd-color-3-s), var(--cd-color-3-l), 0.8)
      }
      .slidecontainer {
        width: 100%;
      }
      .colourSelect {
        background: #eee;
        border-radius: 5px;
        padding: 4px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
      }
      .slider {
        -webkit-appearance: none;
        width: 100%;
        height: 15px;
        background: #d3d3d3;
        border-radius: 5px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
        -webkit-transition: .2s;
        transition: opacity .2s;
      }
      .slider:hover {
        opacity: 1; 
      }
      .slider::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 25px;
        height: 25px;
        border-radius: 50%; 
        background: ${themeColor};
        stroke: dimgrey;
        cursor: pointer;
      }
      .slider::-moz-range-thumb {
        width: 25px;
        height: 25px;
        border-radius: 50%;
        stroke: dimgrey;
        background: ${themeColor};
        cursor: pointer;
      } 
      .tree-container{
        max-height: 1000px;
        overflow: scroll;
      }
      .label{
        display: none;
      }
      .label.show{
        display: inline;
      }
      .node.hovered .label {
          display:inline;
        }
      div.sticky {
          position: -webkit-sticky; /* Safari */
          position: sticky;
          top: 0;
        }
        .searchbar {
          border-style:solid; 
          border-color: lightgrey; 
          border-radius: 5px; 
          float:right
        }
      @media print {
        .tree-container{
        max-height: none;
        overflow: visible;
        
        }
        
        .slider-block {
          display: none;
        }
        .container {
        padding-right: 1.5cm;
        padding-left: 1.5cm;
        padding-bottom: 1.5cm;
        margin: 1cm;
        min-width: 2200px;
        font-size:2.5vw;
        }
        .searchbar {
          display: none;
        }
        h3{ 
          font-size: 2.5vw;
        }
        h2 {
          font-size: 4vw;
          padding: 1cm;
        }
        h1 {
          font-size: 5vw;
        }
        .command-block {
          display: none;
        }
        pre {
          display: none;
        }
        .civet-logo {
          width: 2cm;
          height: 2cm;
        }
        .tree_svg {
          width: 1200px
        }
        .page-footer {
          display: none;
        }
        .civet-header {
          text-align: left;
        }
        .content-block, p {
        page-break-inside: avoid;
        }
      }
      @media screen and (prefers-color-scheme: dark) {
        body {
            background-color: #17141F;
            color: #F2E7DC;
            opacity: 0.95;
          }
          img {
            filter: brightness(.8) contrast(1.2);
          }
        .component {
          background-color: #2C2640;
        }
        .table-striped>tbody>tr:nth-child(odd) {
          background-color: #2C2640;
          /* border-top-color: #3B325B; */
          /* border-color: #2C2640; */
          opacity: 0.9;
        }
        .table {
          border-top: 0px;
          /* border-color: #3B325B; */
          background-color: #17141F;
      }
      .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
        border-top: none;
      }
      .accordion {
          background-color: #2C2640;
          color: #F2E7DC;
          cursor: pointer;
          padding: 12px;
          width: 100%;
          border: none;
          text-align: left;
          outline: none;
          transition: 0.4s;
        }

        .active, .accordion:hover {
          background-color: #17141F;
        }

        .accordion:after {
          content: '\002B';
          color: #F2E7DC;
          font-weight: bold;
          float: right;
          margin-left: 5px;
        }

        .active:after {
          content: "\2212";
        }

        .panel {
          padding: 0 12px;
          background-color: #2C2640;
          max-height: 0;
          overflow: hidden;
          transition: max-height 0.2s ease-out;
        }
      pre {
        background-color: #3B325B;
        color: #F2E7DC;
        border: none;
        opacity: 0.8;
      }
      .searchbar {
        background-color: #3B325B;
        color: #F2E7DC;
        border-style:none;
        opacity: 0.8;
        }
      .slider {
        background: #F2E7DC;
        stroke: #F2E7DC;
      }
      .slider::-webkit-slider-thumb {
        background: #5F9C82;
        fill: #5F9C82;
        stroke: #F2E7DC;
      }
      .slider::-moz-range-thumb {
        stroke: #F2E7DC;
        background: #5F9C82;
        fill: #5F9C82;
      } 
      .node-background{
          fill:#F2E7DC;
          stroke:#F2E7DC;
          opacity: 0.85;
      }
      .node circle{
        stroke-width:0;
        cursor:pointer;
        fill:#5F9C82;
        stroke:#F2E7DC;
        
        }
      .node circle.selected{
        stroke-width:0;
        cursor:pointer;
        fill:#E27E7E;
        stroke:#F2E7DC;
        opacity: 1;
        }
      .node rect{
        stroke-width:2;
        fill:#E27E7E;
        stroke:#F2E7DC;
      }
      .svg-tooltip {
          background: rgba(69,77,93,.9);
          color: #F2E7DC;
    }
    .branch path{
      stroke: #F2E7DC;
      opacity: 0.85;
      }
      .branch.hovered path{
        stroke:#F2E7DC;
        opacity: 1;
      }
        .node.hovered circle{
          stroke:#F2E7DC;
        opacity: 1;
        }
        .node text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
          font-weight: 300;
          font-size: 0.9em;
          color: #F2E7DC;
          fill: #F2E7DC;        
        }
        .scale-bar line {
          stroke: #F2E7DC;
        }
        .scale-bar text{
          fill: #F2E7DC;
          color:  #F2E7DC;
        }
      }
    </style>

  </head>

  <body>
    <script>
      $(document).ready(function() {
        $(window).scroll(function() {
        if ($(this).scrollTop() > 20) {
        $('#toTopBtn').fadeIn();
        } else {
        $('#toTopBtn').fadeOut();
        }
        });
        
        $('#toTopBtn').click(function() {
        $("html, body").animate({
        scrollTop: 0
        }, 400);
        return false;
        });
        });
    </script>
    <!-- <script>
      var colorWell;
      var defaultColor = "#557b86";
      window.addEventListener("load", startup, false);
      function startup() {
        colorWell = document.querySelector("#colorWell");
        colorWell.value = defaultColor;
        colorWell.addEventListener("input", updateFirst, false);
        colorWell.addEventListener("change", updateAll, false);
        colorWell.select();
      }
      function updateFirst(event) {
        var p = document.querySelector("accordion active");
        // var toTopButton = document.getElementById("toTopBtn");
        if (p) {
          p.style.color = event.target.value;
        }
      }
      function updateAll(event) {
        document.querySelectorAll("accordion active").forEach(function(p) {
          p.style.color = event.target.value;
        });
      }
    </script> -->

    <!--Figtree.js-->
    <script type="text/javascript"> 
      const updateTableFactory = (tooltipId,metadata)=>(tipId)=>{
              const data = metadata[tipId];
              const tableDiv = d3.select(document.getElementById(tooltipId));
              //Remove table
          tableDiv.html("")
              if (data !== undefined) {
                  const visibleData = Object.keys(data).filter(d=>d!=='${config["input_display_column"]}');
                  tableDiv.append("h3")
                      .attr("class",'tooltip-header')
                      .text(tipId)
                      .append("hr");
                  tableDiv.selectAll("p")
                          .data(visibleData)
                          .enter()
                          .append("p")
                          .attr("class","tooltip-text")
                              .selectAll("span")
                              .data(d=>[d,data[d]])
                              .enter()
                              .append("span")
                              .attr("class",(d,i)=> i===0? "tooltip-key" :"tooltip-value")
                              .text((d,i)=>i===0? d + " : ": d);
              }
      }
      <%text>
      function addColourEventHandler(circleNodes,legend,colourSelectID,colorCodes,fig){
        d3.select(`#${colourSelectID}`).on("change", function(d){
            const selectedGroup = this.value 
            const colorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations[selectedGroup].values)
            circleNodes.attr("fill",n=>colorScale(n.annotations[selectedGroup]))

            legend.scale(colorScale)
            fig.update();
            console.log(selectedGroup);
          })
        }

      function addTraitColorEventHandler(traits,traitLegend,barSelectID,colorCodes,fig){
        
        d3.select(`#${barSelectID}`).on("change", function(d){
            const selectedGroup = this.value 
            const traitColorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations[selectedGroup].values)
            traits.attr("fill",n => traitColorScale(n.annotations[selectedGroup]))
            traitLegend.scale(traitColorScale)
            fig.update();
            console.log(selectedGroup);
          })
        }


      function addSliderEventHandler(sliderID, fig) {
          const svg = fig.svgSelection.select(function () {
              return this.parentNode;
          })
          const initialHeight = svg.attr("height");
          const maxHeight = fig.tree().externalNodes.length * 50; // 50 pixels for each tip plus a little for margins;
          if (maxHeight <= initialHeight) {
              console.log(sliderID);
              d3.select(`#${sliderID}`).select(function () {
                  return this.parentNode;
              })
                  .remove();// don't need  a slider add names
              fig.svgSelection.selectAll(".label")
                  .classed("show", true)
              return;
          }
          const heightScale = d3.scaleLinear()
              .range([initialHeight, maxHeight])
              .domain([0, 1])
          if (initialHeight / fig.tree().externalNodes.length > 12) {
              fig.svgSelection.selectAll(".label")
                  .classed("show", true)
          }
          d3.select(`#${sliderID}`).on("input", function () {
              const svgHeight = heightScale(this.value);
              //magic number!!
              svg.attr("height", svgHeight);
              fig.update();
              if (svgHeight / fig.tree().externalNodes.filter(node => !node[fig.id].ignore).length > 12) {
                  fig.svgSelection.selectAll(".label")
                      .classed("show", true)
              } else {
                  fig.svgSelection.selectAll(".label")
                      .classed("show", false)
              }
          })
      }
      </%text>
      
      function buildTree(svgID, myTreeString,tooltipID,backgroundDataString,sliderID,colourSelectID,barSelectID,colorCodes) {
          const backgroundData = JSON.parse(backgroundDataString);
          const updateTable = updateTableFactory(tooltipID, backgroundData);
          const margins = {top:50,bottom:60,left:100,right:250}
          const svg = d3.select(document.getElementById(svgID))
          svg.selectAll("g").remove();
          const nexusString = myTreeString;
          const tree = figtree.Tree.parseNexus(nexusString)[0];
          const fig = new figtree.FigTree(document.getElementById(svgID),margins, tree)
          const colorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations["query_boolean"].values)
          const traitColorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations["query_boolean"].values)
          const circleNodes = figtree.circle()
                              .filter(n => !n.children)
                              .attr("r", 8)
                              .attr("fill", n => colorScale(n.annotations["query_boolean"]))
                              .hilightOnHover(20)
                              .onClick((node, i, n) => {
                                  const isSelected = d3.select(n[i]).classed("selected");
                                  fig.svgSelection.selectAll(".selected").classed("selected", false);
                                  if(isSelected){
                                      d3.select(n[i]).classed("selected", false);
                                      updateTable(null);
                                  }else{
                                      d3.select(n[i]).classed("selected", true);
                                      updateTable(node.name);
                                  }
                              });
          const legend = figtree.legend()
                                .scale(colorScale)
                                .x(-100)
                                .y(40)
          const traitLegend = figtree.legend()
                                .scale(traitColorScale)
                                .x(-150)
                                .y(40)
          const traits = figtree.traitBar()
                                .x(svg.style("width")+230)
                                .width(10)
                                .attr("fill",n => traitColorScale(n.annotations["query_boolean"]));
          fig.layout(figtree.rectangularLayout)
                  .nodes(circleNodes,
                          figtree.tipLabel(v=>v.name).attr("dx",10),
                          figtree.rectangle()
                                  .filter(n=>n[fig.id].collapsed)
                                  .attr("width",20)
                                  .attr("height",20)
                  )
                        .nodeBackgrounds(figtree.circle()
                                          .attr("r", 10)
                                .filter(n=>!n.children)
                                        )
                        .branches(figtree.branch()
                                    .hilightOnHover(20) 
                                    .collapseOnClick()
                                    .on("click",()=>{
                                      const svgHeight = fig.svgSelection.select(function() { return this.parentNode; }).attr("height");
                                      if(svgHeight/fig.tree().externalNodes.filter(node=>!node[fig.id].ignore).length>12){
                                        fig.svgSelection.selectAll(".label")
                                          .classed("show",true)
                                      }else{
                                        fig.svgSelection.selectAll(".label")
                                        .classed("show",false)
                                      }
                                    })
                            )
                            .feature(
                                    figtree.scaleBar()
                                      .direction("x")
                                      .length(1/29903)
                                      .x(-60)
                                      .y(-30)
                                      // .y(fig.svgSelection.select(function() { return this.parentNode; }).attr("height")-margins.top-margins.bottom+20)
                                      .title({text:"~1 SNP",
                                      yPadding:10})
                                        )
                            .feature(legend)
                            // .feature(traits)
        addSliderEventHandler(sliderID,fig);
        addColourEventHandler(circleNodes,legend,colourSelectID,colorCodes,fig);
        // addTraitColorEventHandler(traits,traitLegend,barSelectID,colorCodes,fig)
      }
    </script>

<script>
  function exportImageSVG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadSvg(document.querySelector(svgID), name);
      };
  };
  function exportImagePNG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadPng(document.querySelector(svgID), name);
      };
  };
</script>
    <div class="container">
      <a href="#" id="toTopBtn" class="cd-top text-replace js-cd-top cd-top--is-visible cd-top--fade-out" data-abc="true"></a>
      <div>
        <header class="civet-header">
            civet ${version} | 
            <small class="text-muted">Cluster Investigation and Virus Epidemiology Tool</small>
            <hr>
        </header>
        
        <h1>${config["report_title"]}
            <small class="text-muted" style="color:${themeColor}">${date}</small>
        </h1> 
        <br>
        </div>
    <br>
<!--     
    </div>
    <button class="accordion">Report options</button>
    <div class="panel">
      <div class="row">
        <div class="column">
            <input type="color" value="#557b86" id="colorWell">
            <label for="colorWell">Theme Colour</label>
          </div>
        <div class="column">
          </div>
      </div>
    </div> -->
   %if '1' in config['report_content']:
          <h3><strong>Table 1</strong> | Summary of queries </h3>
          <button class="accordion">Table options</button>
          <div class="panel">
            <div class="row">
              <div class="col-sm-2">
                <strong>Show columns:</strong>
              </div>

              <% col_no=0 %>
              %for col in config["query_table_content"]:
                
                <div class="col-sm-1">
                  <a class="toggle-vis" data-column="${col_no}" style="color:${themeColor}">${col.title().replace("_"," ")}</a> 
                </div>
                <% col_no +=1 %>
              %endfor

          </div>
          <div class="row">
            <div class="col-sm-2" ><strong>Export table: </strong></div>
            <div class="col-sm-8" id="tableExportID"></div>
          </div>
          </div>
          <table class="display nowrap" id="myTable">
            <thead>
              <tr>
              %for col in config["query_table_content"]:
              <th style="width:10%;">${col.title().replace("_"," ")}</th>
              %endfor
              </tr>
            </thead>
            <tbody>
              % for row in query_summary_data:
                  <tr>
                    %for col in config["query_table_content"]:
                      %if col=="catchment":
                      <td><a href="#header_${row[col]}" style="color:${themeColor}">${row[col]}</a></td>
                      %else:
                      <td>${row[col]}</td>
                      %endif
                    %endfor
                  </tr>
              % endfor
              </tbody>
            </table>
            
            <script type="text/javascript">
              $(document).ready( function () {
                  var table = $('#myTable').DataTable({
                    "scrollY": "300px",
                    "paging": false,
                    "border-bottom":false,
                    dom: 'frtip',
                    buttons: ["copy","csv","print"]
                  });
                  table.buttons().container().appendTo( $('#tableExportID') );
                  $('a.toggle-vis').on( 'click', function (e) {
                      e.preventDefault();
              
                      // Get the column API object
                      var column = table.column( $(this).attr('data-column') );
              
                      // Toggle the visibility
                      column.visible( ! column.visible() );
                  } );
    
                } );
            </script>
        %endif
        %if 'query_fasta' in config:
    
          <h3><strong>Table 2 </strong> | Queries provided in fasta file</h3>
          <button class="accordion">Passed QC</button>
          <div class="panel">
            <table class="table table-striped" id="myTable2">
            <tr class="header">
                %for col in config["fasta_table_content"]:
                <th style="width:10%;">${col.title().replace("_"," ")}</th>
                %endfor
                </tr>
                % for row in fasta_summary_pass:
                    <tr>
                      %for col in config['fasta_table_content']:
                      <td>${row[col]}</td>
                      %endfor
                    </tr>
                % endfor
            </table>
          </div>
          <button class="accordion">Failed QC</button>
          <div class="panel">
          <table class="table table-striped" id="myTable2">
            <tr class="header">
                %for col in config["fasta_table_content"]:
                <th style="width:10%;">${col.title().replace("_"," ")}</th>
                %endfor
                </tr>
                % for row in fasta_summary_fail:
                    <tr>
                      %for col in config['fasta_table_content']:
                      <td>${row[col]}</td>
                      %endfor
                    </tr>
                % endfor
            </table>
          </div>
        %endif
    <% figure_count = 0 %>
        %if config["global_snipit"]:
          <h2>global snipit </h2>
          <h3>SNP summary for all query sequences relative to the reference </h3>
          <% figure_count +=1 %>
                <br>
                <button class="accordion">Export image</button>
                  <div class="panel">
                    <div class="row">
                      <div class="column">
                        <button id="global_snipit_svg">SVG</button>
                      </div>
                      <div class="column">
                        <button id="global_snipit_png">PNG</button>
                      </div>
                    </div>
                  </div>
                    <div id="global_snipit">
                    ${data_for_report["global_snipit_svg"]}
                    </div>
              <script type="text/javascript">
                exportImageSVG("#global_snipit_svg","#global_snipit", "global_snipit_chart");
              </script>
              <script type="text/javascript">
                exportImagePNG("#global_snipit_png","#global_snipit", "global_snipit_chart");
              </script>
                    <h3><strong>Figure ${figure_count}</strong> | snipit plot for all focal query sequences</h3>
                    <hr>        
          %endif
        %if '8' in config["report_content"]:
        <%timeseries_data = data_for_report["full_metadata"] %>
        <%date_field = config["input_date_column"] %>
        <%series_colour_factor = config["series_colour_factor"] %>

        <% figure_count +=1 %>
        <br>

        <div id="time_series" style="width:95%"></div>

          <script>
            var vlSpec_time = {
              "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
              "width": "container",
              "height": 200,
              "datasets": {"time_series": ${timeseries_data}},
              "data": {
                "name": "time_series"
                  },
                "mark": "bar",
                "encoding": {
                  "x": {
                    "field": "${date_field}", 
                    "bin":true,
                    "scale": {"type": "nice","interval": "week", "step": 2},
                    "timeUnit": {
                        "unit": "utcyearmonthdate",
                        "step": 7},
                    "title": "Date",
                    "axis": {
                      "grid": false,
                      "format":"%Y-%m-%d",
                      "labelFont":"Helvetica Neue",
                      "labelFontSize":18,
                      "titleFontSize":18,
                      "titleFont":"Helvetica Neue"
                    },
                  },
                  "y": 
                  {"aggregate": "count",
                  "title": "Genome count",
                  "axis":{
                  "grid": false,
                  "labelFont":"Helvetica Neue",
                  "labelFontSize":18,
                  "titleFontSize":18,
                  "titleFont":"Helvetica Neue"}
                  },
                  "color": {
                        "field": "${series_colour_factor}", 
                        "type": "nominal",
                        "scale": {
                              "range": [
                                    "#B6B8C8",
                                    "#D4B489",
                                    "#A6626F",
                                    "#733646",
                                    "#A47E3E",
                                    "#DC9598",
                                    "#83818F",
                                    "#B3ABD0",
                                    "#B8B2C4",
                                    "#A07E62",
                                    "#F9C0C7"
                                  ]
                            },
                        "legend": {
                            "title": "${series_colour_factor.capitalize()}",
                            "labelFontSize":14,
                            "labelFont":"Helvetica Neue",
                            "titleFontSize":16,
                            "titleFont":"Helvetica Neue",
                            "titleFontStyle":"normal"
                            }
                        }
                    },
                        "config": {
                          "view": {"stroke": null},
                          "axis": {"grid": false},
                          "text": {"font":"Helvetica Neue","fontWeight":0.1}
                        }
                };          
          vegaEmbed('#time_series', vlSpec_time, {renderer: "svg"})
                .then(result => console.log(result))
                .catch(console.warn);
  </script>
          <% figure_count +=1 %>
          <h3><strong>Figure ${figure_count}</strong> | Time series of query sequences and epidata</h3>
          
          <div id="time_series_fig" style="width:90%"></div>
  %endif

    %for catchment in catchments:
        <% catchment_name = catchment.replace("_"," ").title() %>
        <h2><a id = "header_${catchment}"></a>${catchment_name}</h2> 
        
        %if '2' in config["report_content"]:
          %if 'query_fasta' in config:
            <h3><strong>Table ${2+int(catchment_name.split(" ")[1])}</strong> | Summary of ${catchment_name}</h3>
          %else:
            <h3><strong>Table ${1+int(catchment_name.split(" ")[1])}</strong> | Summary of ${catchment_name}</h3>
          %endif
            <table class="table table-striped" id='${data_for_report[catchment]["catchmentID"]}'>
                  <tr class="header">
                  </tr>
                  <table class="table table-striped">
                  <tr class="header">
                  <th style="width:30%;">Information</th>
                  <th style="width:60%;"> </th>
                  </tr>
                  % for stat in data_for_report[catchment]['catchment_summary_data']:
                    <tr>
                      <td>${stat.title().replace("_"," ")}</td>
                      <td>${data_for_report[catchment]['catchment_summary_data'][stat]}</td>
                    </tr>
                  % endfor
                </table>
            %endif

        %if '3' in config["report_content"] and catchment in config["figure_catchments"]:
         <% figure_count +=1 %>
            <button class="accordion">Tree options</button>
            <div class="panel">
              <div class="row">
                <div class="column">
                  <div class="slider-block" id="slider_${data_for_report[catchment]['catchmentID']}">
                    <p>Expansion</p>
                    <input class="slider" type="range" id="rangeinput_${data_for_report[catchment]['catchmentID']}"  min="0" max="1" style="width: 100px" step="0.01" value="0" />
                    <span class="highlight"></span>
                  </div>
                </div>
                <div class="column">
                  <div>
                  <p>Colour by</p>
                  <select class="colourSelect" id="colourSelect_${data_for_report[catchment]['catchmentID']}">
                    <option value="query_boolean">Query</option>
                    % for annotation in config['tree_annotations'].split(" "):
                      <option value="${annotation}">${annotation.title()}</option>
                    % endfor
                  </select>
                  </div>
                </div>
                <!-- <div class="column">
                  <div>
                  <p>Colour by</p>
                  <select class="colourSelect" id="barSelect_${data_for_report[catchment]['catchmentID']}">
                    <option value="query_boolean">Query</option>
                    <option value="country">Country</option>
                    <option value="lineage">Lineage</option>
                  </select>
                  </div>
                </div> -->
              </div>
            </div>
      
        
          <div class="row tree-container">

            <div class="col-xs-7">
              <svg class="tree_svg" width="700" height="400" id="tree_${data_for_report[catchment]['catchmentID']}"></svg>
            </div>
            <div class="col-xs-4 sticky" id="tooltip_${data_for_report[catchment]['catchmentID']}">
            </div> 
            
            <script type="text/javascript">
              buildTree("tree_${data_for_report[catchment]['catchmentID']}", 
                        `${data_for_report[catchment]['nexus']}`,
                        "tooltip_${data_for_report[catchment]['catchmentID']}",
                        '${background_data}',
                        "rangeinput_${data_for_report[catchment]['catchmentID']}",
                        "colourSelect_${data_for_report[catchment]['catchmentID']}",
                        "barSelect_${data_for_report[catchment]['catchmentID']}",
                        ${colorCodes});
            </script> 
          </div> 
          <h3><strong>Figure ${figure_count}</strong> | ${catchment_name} phylogeny</h3>
          <hr>
        
        %endif
        %if '4' in config["report_content"] and catchment in config["figure_catchments"]:
        <% figure_count +=1 %>
        <br>
        <button class="accordion">Export image</button>
          <div class="panel">
            <div class="row">
              <div class="column">
                <button id="${catchment}_snipit_svg">SVG</button>
              </div>
              <div class="column">
                <button id="${catchment}_snipit_png">PNG</button>
              </div>
            </div>
          </div>
            <div id="${catchment}_snipit">
            ${data_for_report[catchment]["snipit_svg"]}
            </div>
      <script type="text/javascript">
        exportImageSVG("#${catchment}_snipit_svg","#${catchment}_snipit","${catchment}_snipit_graph");
      </script>
      <script type="text/javascript">
        exportImagePNG("#${catchment}_snipit_png","#${catchment}_snipit","${catchment}_snipit_graph");
      </script>
            <h3><strong>Figure ${figure_count}</strong> | snipit plot for queries in ${catchment_name}</h3>
            <hr>
            
        %endif
        %if '5' in config["report_content"] and catchment in config["figure_catchments"]:
        <%timeline_data = data_for_report[catchment]["timeline_data"] %>
        <% timeline_height = 50+(30 * data_for_report[catchment]['catchment_summary_data']["query_count"]) %>
        
        <% figure_count +=1 %>
        <br>
              <div id="${catchment}_timeline" style="width:90%"></div>
                <script>
                  var vlSpec_time = {
                    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
                    "width": "container",
                    "height": ${timeline_height},
                    "datasets": ${timeline_data},
                    "data": {
                      "name": "${catchment}_timeline"
                        },
                      "transform":[
                        {
                        "filter": {
                          "field": "date_type",
                          "oneOf": ${config["timeline_dates"]}
                        }
                        }
                      ],
                      "encoding": {
                        "x": {
                          "field": "date", 
                          "timeUnit":"monthdate",
                          "type": "temporal",
                          "title": "Date",
                          "axis": {
                            "grid": false,
                            "labelFont":"Helvetica Neue",
                            "labelFontSize":18,
                            "tickCount": {"interval": "day", "step": 2}, // put in num days/
                            "titleFontSize":18,
                            "titleFont":"Helvetica Neue"
                          },
                        },
                        "y": 
                        {"field": "sequence_name", 
                        "type": "nominal",
                        "axis": {
                          "title": "",
                          "font": "Helvetica Neue",
                          "labelFontSize":15,
                          "labelLimit": 0,
                          "labelPadding":20
                        }
                      },
                       "color": {
                          "field": "date_type", 
                          "type": "nominal",
                          "legend": {
                            "labelFont":"Helvetica Neue",
                            "labelFontSize":18,
                            "orient": "bottom", 
                            "symbolSize": 500,
                            "title": null}
                        }
                      },
                        "config": {
                          "view": {"stroke": null},
                          "axis": {"grid": false},
                          "text": {"font":"Helvetica Neue"}
                        },
                      "layer": [
                        {
                        "mark": "line",
                        "encoding": {
                          "detail": {
                          "field": "sequence_name",
                          "type": "nominal"
                          },
                          "color": {"value": "#A9A9A9"}
                        }
                        },
                        {
                        "mark": {
                          "type": "point",
                          "filled": true,
                          "tooltip":true
                        },
                        "encoding": {
                          "tooltip":[
                            {"field":"date","type":"temporal","title":"Date"}
                        ],
                          "color": {
                          "field": "date_type",
                          "type": "nominal",
                          "scale": {
                            "domain": ${config["timeline_dates"]},
                            "range": ${colorCodes}
                          },
                          "title": "Date"
                          },
                          "size": {"value": 500},
                          "opacity": {"value": 1}
                        }
                        }
                      ]
                      };          
                vegaEmbed('#${catchment}_timeline', vlSpec_time, {renderer: "svg"})
                      .then(result => console.log(result))
                      .catch(console.warn);
              </script>
        <h3><strong>Figure ${figure_count}</strong> | Timeline plot for queries in ${catchment_name}</h3>
        <hr>
        %endif
        %endfor
        %if '6' in config["report_content"]:
          <% figure_count +=1 %>
          <br>
        <div id="background_map" style="width:90%"></div>
        <script>
        var vSpec_bmap = {
          "$schema": "https://vega.github.io/schema/vega/v5.json",
          "autosize": {"type":"fit", "contains":"padding"},
          "background": "white",
          "padding": 5,
          "width": 500,
          "height": 500,
          "style": "cell",

% if config["civet_mode"] == "CLIMB":
    <%start_scale = 750%>
    <%max_scale = 10000%>
    <%start_x = 4%>
    <%start_y = 55%>
% else:
    <%start_scale = 150%>
    <%max_scale = 4000%>
    <%start_x = 0%>
    <%start_y = 0%>
% endif

          "signals": [
                      { "name": "tx", "update": "width / 2" },
                      { "name": "ty", "update": "height / 2" },
                      {
                        "name": "scale",
                        "value": ${start_scale},
                        "on": [{
                          "events": {"type": "wheel", "consume": true},
                          "update": "clamp(scale * pow(1.0005, -event.deltaY * pow(16, event.deltaMode)), 100, ${max_scale})"
                        }]
                      },
                      {
                        "name": "angles",
                        "value": [0, 0],
                        "on": [{
                          "events": "mousedown",
                          "update": "[rotateX, centerY]"
                        }]
                      },
                      {
                        "name": "cloned",
                        "value": null,
                        "on": [{
                          "events": "mousedown",
                          "update": "copy('projection')"
                        }]
                      },
                      {
                        "name": "start",
                        "value": null,
                        "on": [{
                          "events": "mousedown",
                          "update": "invert(cloned, xy())"
                        }]
                      },
                      {
                        "name": "drag", "value": null,
                        "on": [{
                          "events": "[mousedown, window:mouseup] > window:mousemove",
                          "update": "invert(cloned, xy())"
                        }]
                      },
                      {
                        "name": "delta", "value": null,
                        "on": [{
                          "events": {"signal": "drag"},
                          "update": "[drag[0] - start[0], start[1] - drag[1]]"
                        }]
                      },
                      {
                        "name": "rotateX", "value": ${start_x},
                        "on": [{
                          "events": {"signal": "delta"},
                          "update": "angles[0] + delta[0]"
                        }]
                      },
                      {
                        "name": "centerY", "value": ${start_y},
                        "on": [{
                          "events": {"signal": "delta"},
                          "update": "clamp(angles[1] + delta[1], -60, 60)"
                        }]
                      }
%for location in data_for_report['locations_wanted']:
<%
start_arc = data_for_report[location]["start_arc"]
start_inner_arc = data_for_report[location]["start_inner_arc"]
start_text = data_for_report[location]["start_text"]%>
                      ,{
                        "name": "arc_zoom_${location}",
                        "value": ${start_arc}
                        },{
                        "name":"inner_arc_zoom_${location}",
                        "value":${start_inner_arc}
                        
                    },{
                      "name":"text_zoom_${location}",
                      "value":${start_text}
                    }
  %endfor
                    ],
          "data": [
            {
              "name": "background_data",
              "url": "${config['background_map_file']}",
              "format": {"type": "topojson", "feature": "${config['background_topojson_feature_name']}"}
            },
            {
              "name": "lineage_data",
              "values": 
                ${data_for_report["background_map_data"]}
              ,
              "format": {}
            }
%for location in data_for_report['locations_wanted']:
<%latitude = data_for_report[location]["centroids"][1]
longitude = data_for_report[location]["centroids"][0]%>
            ,{
              "name": "data_${location}_1",
              "source": "lineage_data",
              "transform": [{"type": "filter", "expr": "datum.location == '${location}'"}]
            },
            {
              "name": "data_${location}_2",
              "source": "data_${location}_1",
              "transform": [
                {
                  "type": "geopoint",
                  "projection": "projection",
                  "fields": [{"expr": "${longitude}"}, {"expr": "${latitude}"}],
                  "as": ["layer_${location}_layer_0_x", "layer_${location}_layer_0_y"]
                },
                {
                  "type": "stack",
                  "groupby": [],
                  "field": "count",
                  "sort": {"field": ["count"], "order": ["ascending"]},
                  "as": ["count_start", "count_end"],
                  "offset": "zero"
                },
                {
                  "type": "filter",
                  "expr": "isValid(datum[\"count\"]) && isFinite(+datum[\"count\"])"
                }
              ]
            },
            {
              "name": "data_${location}_3",
              "source": "data_${location}_1",
              "transform": [
                {
                  "type": "geopoint",
                  "projection": "projection",
                  "fields": [{"expr": "${longitude}"}, {"expr": "${latitude}"}],
                  "as": ["layer_${location}_layer_1_x", "layer_${location}_layer_1_y"]
                },
                {
                  "type": "stack",
                  "groupby": [],
                  "field": "count",
                  "sort": {
                    "field": ["count", "count"],
                    "order": ["ascending", "ascending"]
                  },
                  "as": ["count_start", "count_end"],
                  "offset": "zero"
                },
                {
                  "type": "filter",
                  "expr": "isValid(datum[\"count\"]) && isFinite(+datum[\"count\"])"
                }
              ]
            }
%endfor
          ],

          "projections": [
            {
            "name": "projection",
            "type": "mercator",
            "scale": {"signal": "scale"},
            "rotate": [{"signal": "rotateX"}, 0, 0],
            "center": [0, {"signal": "centerY"}],
            "translate": [{"signal": "tx"}, {"signal": "ty"}]
            }
          ],
          "marks": [
            {
              "name": "layer_0_marks",
              "type": "shape",
              "clip": true,
              "style": ["geoshape"],
              "from": {"data": "background_data"},
              "encode": {
                "enter": {
                    "strokeWidth": {"value": 0.2},
                    "stroke": {"value": "white"},
                    "fill": {"value": "lightgrey"}
                  }
              },
              "transform": [{"type": "geoshape", "projection": "projection"}]
            }
%for location in data_for_report['locations_wanted']:
            ,{
              "name": "layer_${location}_layer_0_marks",
              "type": "arc",
              "clip": true,
              "style": ["arc"],
              "from": {"data": "data_${location}_2"},
              "encode": {
                "update": {
                  "stroke": {"value": "white"},
                  "strokeWidth": {"value": 0.5},
                  "padAngle": {"value": 0.05},
                  "innerRadius": {"signal": "inner_arc_zoom_${location}"},
                  "tooltip": {
                    "signal": "{\"location\": datum[\"location\"], \"lineage\": isValid(datum[\"lineage\"]) ? datum[\"lineage\"] : \"\"+datum[\"lineage\"], \"count\": format(datum[\"count\"], \"\")}"
                  },
                  "cornerRadius": {"value": 5},
                  "fill": {"scale": "color", "field": "lineage"},
                  "description": {
                    "signal": "\"count: \" + (format(datum[\"count\"], \"\")) + \"; lineage: \" + (isValid(datum[\"lineage\"]) ? datum[\"lineage\"] : \"\"+datum[\"lineage\"])"
                  },
                  "x": {"field": "layer_${location}_layer_0_x"},
                  "y": {"field": "layer_${location}_layer_0_y"},
                  "outerRadius": {"scale": "radius_${location}", "field": "count"},
                  "startAngle": {
                    "scale": "theta_${location}",
                    "field": "count_end",
                    "offset": -1.5
                  },
                  "endAngle": {"scale": "theta_${location}", "field": "count_start", "offset": -1.5}
                }
              }
            },
            {
              "name": "layer_${location}_layer_1_marks",
              "type": "text",
              "clip": true,
              "style": ["text"],
              "from": {"data": "data_${location}_3"},
              "encode": {
                "update": {
                  "fill": {"value": "black"},
                  "stroke": {"value": "white"},
                  "strokeWidth": {"value": 0.1},
                  "align": {"value": "center"},
                  "radius": {"scale": "radius_${location}", "field": "count","offset":10},
                  "font": {"value": "Helvetica Neue"},
                  "fontSize": {"signal": "text_zoom_${location}"},
                  "description": {
                    "signal":  "\"data: \" + (format(datum[\"data\"], \"\"))"},
                  "x": {"field": "layer_${location}_layer_1_x"},
                  "y": {"field": "layer_${location}_layer_1_y"},
                  "text": {
                    "signal": "isValid(datum[\"lineage\"]) ? datum[\"lineage\"] : \"\"+datum[\"lineage\"]"
                  },
                  "baseline": {"value": "middle"},
                  "theta": {
                    "signal": "scale(\"theta_${location}\", 0.5 * datum[\"count_start\"] + 0.5 * datum[\"count_end\"])",
                    "offset": -1.5
                  }
                }
              },
              "transform": [
                {"type": "label", "size": [800, 500]}
              ]
              }
%endfor
          ]
          ,
          "scales": [
%for location in data_for_report['locations_wanted']:
            {
              "name": "theta_${location}",
              "type": "linear",
              "domain": {
                "fields": [
                  {"data": "data_${location}_2", "field": "count_start"},
                  {"data": "data_${location}_2", "field": "count_end"},
                  {"data": "data_${location}_3", "field": "count_start"},
                  {"data": "data_${location}_3", "field": "count_end"}
                ]
              },
              "range": [0, 6.283185307179586],
              "reverse": true,
              "zero": true
            },
            {
              "name": "radius_${location}",
              "type": "sqrt",
              "domain": {
                "fields": [
                  {"data": "data_${location}_2", "field": "count"},
                  {"data": "data_${location}_3", "field": "count"}
                ]
              },
              "range": [{"signal":"inner_arc_zoom_${location}"},{"signal":"arc_zoom_${location}"}],
              "zero": true
            },
%endfor
            {
              "name": "color",
              "type": "ordinal",
              "domain": {
                "fields": [
%for count,location in enumerate(data_for_report['locations_wanted']):
                  {"data": "data_${location}_2", "field": "lineage"},
                  {"data": "data_${location}_3", "field": "lineage"}
%if count < len(data_for_report['locations_wanted'])-1:
,
%endif
%endfor
                ],
                "sort": true
              },
              "range": [
                "#B6B8C8",
                "#D4B489",
                "#A6626F",
                "#733646",
                "#A47E3E",
                "#DC9598",
                "#83818F",
                "#B3ABD0",
                "#B8B2C4",
                "#A07E62",
                "#F9C0C7"
              ]
            }
          ]
        }

      vegaEmbed('#background_map', vSpec_bmap, {renderer: "svg"})
                      .then(result => console.log(result))
                      .catch(console.warn);

      </script>
      <h3><strong>Figure ${figure_count}</strong> | Background diversity map</h3>
      <hr>
        %endif
        
        %if '7' in config["report_content"]:
        <%qmap_data = data_for_report["query_map_data"] %>
        <% figure_count +=1 %>
        <div id="query_map" style="width:90%"></div>
        <script>
			      var vSpec_qmap = {

            "$schema": "https://vega.github.io/schema/vega/v5.json",
            "autosize": "none",
            "background": "white",
            "padding": 5,
            "width": 100,
            "height": 100,

          "signals": [
              { "name": "tx", "update": "width / 2" },
              { "name": "ty", "update": "height / 2" },
              {
                "name": "scale",
                "value": 1000,
                "on": [{
                  "events": {"type": "wheel", "consume": true},
                  "update": "clamp(scale * pow(1.0005, -event.deltaY * pow(16, event.deltaMode)), 100, 3000)"
                }]
              },
              {
                "name": "angles",
                "value": [0, 0],
                "on": [{
                  "events": "mousedown",
                  "update": "[rotateX, centerY]"
                }]
              },
              {
                "name": "cloned",
                "value": null,
                "on": [{
                  "events": "mousedown",
                  "update": "copy('projection')"
                }]
              },
              {
                "name": "start",
                "value": null,
                "on": [{
                  "events": "mousedown",
                  "update": "invert(cloned, xy())"
                }]
              },
              {
                "name": "drag", "value": null,
                "on": [{
                  "events": "[mousedown, window:mouseup] > window:mousemove",
                  "update": "invert(cloned, xy())"
                }]
              },
              {
                "name": "delta", "value": null,
                "on": [{
                  "events": {"signal": "drag"},
                  "update": "[drag[0] - start[0], start[1] - drag[1]]"
                }]
              },
              {
                "name": "rotateX", "value": 0,
                "on": [{
                  "events": {"signal": "delta"},
                  "update": "angles[0] + delta[0]"
                }]
              },
              {
                "name": "centerY", "value": ${config["start_centre_lat"]},
                "on": [{
                  "events": {"signal": "delta"},
                  "update": "clamp(angles[1] + delta[1], -60, 60)"
                }]
              }
            ],

            "projections": [
              {
                "name": "projection",
                "type": "mercator",
                "scale": {"signal": "scale"},
                "rotate": [{"signal": "rotateX"}, 0, 0],
                "center": [${config["start_centre_long"]}, {"signal": "centerY"}],
                "translate": [{"signal": "tx"}, {"signal": "ty"}]
              }
            ],
            "scales":[
              {
                "name":"color",
                "type":"ordinal",
                "domain":{"data":"source_0", "field":"catchment"},
                "range":"category"
              }
            ],

            "data": [
              {
                "name": "background_data",
                "url": "${config['query_map_file']}",
                "format": {"type": "topojson", "feature": "${config['query_topojson_feature_name']}"}
              },
              {
                "name": "source_0",
                "values": ${qmap_data}
              },
              {
                "name": "data_0",
                "source": "background_data"
              },
              {
                "name": "data_1",
                "source": "source_0",
                "transform": [
                  {
                    "type": "geopoint",
                    "projection": "projection",
                    "fields": ["longitude", "latitude"],
                    "as": ["layer_1_x", "layer_1_y"]
                  }
                ]
              }
            ],
            "marks": [
               {
                "name":"layer_0_marks",
                "type": "shape",
                "from": {"data": "background_data"},
                "encode": {
                  "enter": {
                    "strokeWidth": {"value": 0.2},
                    "stroke": {"value": "white"},
                    "fill": {"value": "lightgrey"}
                  }
                },
                "transform": [{"type": "geoshape", "projection": "projection"}]
              },
              {
                "name": "layer_2_marks",
                "type": "symbol",
                "clip": true,
                "style": ["circle"],
                "from": {"data": "data_1"},
                "encode": {
                  "update": {
                    "opacity": {"value": 0.7},
                    "tooltip": {"signal": "datum"},
                    "fill":{"value":"black"},
                    "ariaRoleDescription": {"value": "circle"},
                    "description": {
                      "signal": "\"longitude: \" + (format(datum[\"longitude\"], \"\")) + \"; latitude: \" + (format(datum[\"latitude\"], \"\"))"
                    },
                    "x": {"field": "layer_1_x"},
                    "y": {"field": "layer_1_y"},
                    "size": {"value": 15},
                    "shape": {"value": "circle"}
                  }
                }
              },
              {
                "name": "layer_1_marks",
                "type": "symbol",
                "clip": true,
                "style": ["circle"],
                "from": {"data": "data_1"},
                "encode": {
                  "update": {
                    "opacity": {"value": 0.7},
                    "tooltip": {"signal": "datum"},
                    "fill":{"scale":"color","field":"catchment"},
                    "ariaRoleDescription": {"value": "circle"},
                    "description": {
                      "signal": "\"longitude: \" + (format(datum[\"longitude\"], \"\")) + \"; latitude: \" + (format(datum[\"latitude\"], \"\"))"
                    },
                    "x": {"field": "layer_1_x"},
                    "y": {"field": "layer_1_y"},
                    "size": {"value": 10},
                    "shape": {"value": "circle"}
                  }
                }
              }
            ]
      }
      vegaEmbed('#query_map', vSpec_qmap, {renderer: "svg"})
                      .then(result => console.log(result))
                      .catch(console.warn);
        </script>
      <h3><strong>Figure ${figure_count}</strong> | Query map</h3>
      <hr>
      %endif

    <script>
        var acc = document.getElementsByClassName("accordion");
        var i;
        for (i = 0; i < acc.length; i++) {
              acc[i].addEventListener("click", function() {
                this.classList.toggle("active");
                var panel = this.nextElementSibling;
                if (panel.style.maxHeight) {
                  panel.style.maxHeight = null;
                } else {
                  panel.style.maxHeight = panel.scrollHeight*1.2 + "px";
                } 
              });
            }
    </script>

    <footer class="page-footer">
      <div class="container-fluid text-right text-md-right">
        <hr>
        <div class="row">
          <!--<div class="col-sm-1">
            <p>
            <img class="civet-logo" src="https://raw.githubusercontent.com/COG-UK/civet/master/docs/doc_figures/virus.svg" vertical-align="left" width="50" height="50"></img>
            <p> -->
        </div>

      <div class="col-sm-11" style="text-align: right;">
        civet ${version} | <small class="text-muted">Cluster Investigation and Virus Epidemiology Tool</small> <br><small class="text-muted">GNU General Public License v3.0</small></div>

        <br><br>
        </p>
      </div>
    </footer>
    </div>
  </body>
</html>

function init( json, tabContainer ) {
    d3.json( json, function(err,json) {
        if( err ) { alert("Error loading json" ); }
        d3.select( 'body' )
            .append('div')
            .attr('id',tabContainer)
            .attr('width','100%')
            .attr('height','100%')
        var cg = new CoverageGraph(json,tabContainer);
    });
}

/*
 * Represents the entire page containing tabs for every reference
 * in the given json. Json can be generated using the bam_to_qualdepth.py
 * script.
*/
function CoverageGraph( refjson, tabContainer ) {
    var self = this;
    self.tabs;
    self.tabContainerElement = tabContainer;
    self.tabContainerId = '#'+tabContainer;
    self.json = refjson;

    createTabs( Object.keys(refjson) );

    function createTabs( labels ) {
        //labels.unshift('testtab');

        // Create the tabs
        d3.select(self.tabContainerId)
            .append( 'ul' )
            .attr( 'id', 'tabs' )
            .selectAll('li')
            .data( labels )
            .enter()
            .append( 'li' )
            .attr( 'id', function(d){return d;} )
            .append( 'a' )
            .property( 'href', function(d,i){return '#tabs-' + i;} )
            .style( 'font-size', '8pt' )
            .text( function(d){ return d; } );


        // Create the panels
        d3.select(self.tabContainerId)
            .selectAll( 'div' )
            .data( labels )
            .enter()
            .append( 'div' )
            .attr( 'id', function(d,i){return 'tabs-' + i;} )

        $(self.tabContainerId).tabs({
            create: function( event, ui ) {
                var selectedTabIndex = $(self.tabContainerId).tabs('option','active');
                var selectedTabName = $($(self.tabContainerId+' li')[selectedTabIndex]).text();
                var selectedPanelId = $($(self.tabContainerId+' div')[selectedTabIndex]).attr('id');
                loadTab( selectedTabName, selectedPanelId );
            },
            activate: function( event, ui ) {
                loadTab( ui.newTab.attr('id'), ui.newPanel.attr('id') );
            }
        });
    }

    self.leftAxisScale;
    self.rightAxisScale;
    self.bottomAxisScale;
    function loadTab( refname, panelId ){
        console.log( "Selected Tab Name: " + refname );
        console.log( "Active Panel ID: " + panelId );
        var refinfo = self.json[refname];
        var w = 1366;
        var h = 768;
        var svgId = 'svg-'+refname;
        self.leftAxisScale = d3.scale.linear()
            .domain([0,refinfo['maxd']])
            .range([10, h])
        self.rightAxisScale = d3.scale.linear()
            .domain([0,refinfo['maxq']])
            .range([0, h/3])
        self.bottomAxisScale = d3.scale.linear()
            .domain([0,refinfo['length']])
            .range([0,w])

        // Create the svg
        var svg = d3.select('#'+panelId)
            .append('svg')
            .attr( 'width', w )
            .attr( 'height', h )
            .attr( 'id', svgId )

        /* Left axis position */
        var xleft = self.bottomAxisScale(200);
        /* Width of xmarks */
        var linew = self.bottomAxisScale(100);
        /* Right axis position */
        var xright = self.bottomAxisScale(w-500);

        /* Magic to make the axis */
        makeLeftAxis( svg, xleft, linew, refinfo['maxd'] );
        makeBottomAxis( svg, xleft, linew, refinfo['length'], h );
    }

    function makeBottomAxis( svg, xleft, linew, maxX, maxY ) {
        var maxlefty = self.leftAxisScale(maxY);
        var bottomVals = d3.range(0,maxX,maxX/10);
        bottomVals.push(maxX);
        var bottomTicks = svg.append('g')
            .attr('id','bottomTicks')
            .selectAll('line')
            .data(bottomVals)
            .enter()
            .append('g')
            .each(function(d,i){
                var g = d3.select(this);
                var x = self.bottomAxisScale(d)+xleft+(linew/2);
                g.append('line')
                    .attr('y1',maxY-linew)
                    .attr('y2',maxY-linew*3)
                    .attr('x1',x)
                    .attr('x2',x)
                    .attr('stroke','black')

                g.append('text')
                    .attr('y',maxY)
                    .attr('x',x-10)
                    .attr('font-size',8)
                    .text(d)
            })
            .append('line')
            .attr('x1',xleft)
            .attr('x2',maxX)
            .attr('y1',maxY-linew*2)
            .attr('y2',maxY-linew*2)
            .attr('stroke','black')
    }

    /*
     * Creates the Left Axis with marks and labels
    */
    function makeLeftAxis( svg, xleft, linew, maxY ) {
        var maxlefty = self.leftAxisScale(maxY);
        var leftVals = d3.range(0,maxY,maxY/10);
        leftVals.push( maxY);
        var leftTicks = svg.append('g')
            .attr('id','leftTicks')
            .selectAll('line')
            .data( leftVals )
            .enter()
            .append('g') // Contains the left axis
            .each(function(d,i){
                /* Draw all the individual left x ticks and labels */
                var g = d3.select(this);
                var y = self.leftAxisScale(d)-linew;
                g.append('line')
                    .attr('y1',y-linew)
                    .attr('y2',y-linew)
                    .attr('x1',xleft)
                    .attr('x2',xleft+linew)
                    .attr('stroke','black');

                g.append('text')
                    .attr('y',maxlefty-linew*3-(y-10))
                    .attr('x',0)
                    .attr('font-size',8)
                    .text(d);
            })
            .append( 'line' )
            .attr('x1', xleft+(linew/2))
            .attr('x2', xleft+(linew/2))
            .attr('y1', 0 )
            .attr('y2', maxlefty-linew*2 )
            .attr('stroke','black')
    }
}

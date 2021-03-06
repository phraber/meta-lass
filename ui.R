library(shiny)

logo.formats       <- c("pdf", "eps", "png", "png_print", "svg", "jpeg")
names(logo.formats) = c("pdf", "eps", "png", "hi-res png", "svg", "jpeg")
#if (!system.file("gs") & !system.file("")

shinyUI(

    fluidPage(
        tags$head(
        ),

#        titlePanel("LASSIE: Longitudinal Antigenic Swarm Selection from Intrahost Evolution"),

        tabsetPanel(id="tabset", type="tabs",

            tabPanel("Step 1: Select Sites", id="sites",

                sidebarLayout(

                    sidebarPanel(

                        strong("Input Options"),

                        checkboxInput('demo_data', 'Use CH505 example data',
                            FALSE), #TRUE),

                        c("or choose your own"),

                        fileInput("aas_file", "Protein Alignment File"),

                        sliderInput("tf_loss_cutoff",
                            "Exclude sites with % peak TF loss below:",
                            min = 0, max = 100, value = 100),
hr(),
                        checkboxInput('sites_advanced',
                            'Show more options', FALSE),

                        conditionalPanel("input.sites_advanced",

                        sliderInput("tf_loss_when_up",
                            "Order sites by when they first reach this % TF loss:",
                            min = 0, max = 100, value = 10),

                            radioButtons("aln.format", "Alignment Format",
                                c("fasta", "clustal", "phylip", "msf", "mase"),
                                inline=T )),

                        conditionalPanel("input.sites_advanced",
                            checkboxInput('map.pngs2o',
                                'Map N to O in Nx[ST] PNG motifs.',
                                TRUE)),

                        conditionalPanel("input.sites_advanced",
                            p(strong("Timepoint Parsing"))),

                        conditionalPanel("input.sites_advanced",
                            radioButtons("tp.sep", "Field separator",
                                c(".", "_", "-", "|"), inline=T )),

                        conditionalPanel("input.sites_advanced",
                            numericInput("tp.pos",
                                "Which field contains timepoint",
                                value=1, min=1, max=6, step=1))#,
#hr(),
#                        strong("Review selected sites then check"),
#                        checkboxInput('select_seqs',
#                            'Send Selected Sites to Step 2', FALSE)
                    ),

                    mainPanel(

                        strong("Results: Selected Sites"),

                        p(textOutput("nSelectedSites")),

			dataTableOutput("selectedSites"),

			conditionalPanel("output.selectedSites",

			    hr(),
                            p(strong("Review timepoint parsing"),
			    checkboxInput('by_timepoint',
				'Tabulate sequence counts by timepoint', FALSE)),
			    conditionalPanel("input.by_timepoint",
			        tableOutput("tabulateByTimepoint")),

#			conditionalPanel("output.selectedSites",
                            hr(),
                            p(strong("Review selected sites and check the box below"),
				checkboxInput('select_seqs',
				    'Send selected sites to Step 2', FALSE)
                            ),
                            hr()
                        ) # end of conditional panel
                    ) # end of main panel
                ) # end of sidebar layout
            ), # end of tab panel

            tabPanel("Step 2: Select Sequences", id="seqs",

                sidebarLayout(

                    sidebarPanel(

                        strong(p(textOutput("showSummary"))),

                        conditionalPanel("input.select_seqs",

                            checkboxInput('seqs_advanced',
                                'Show more options', FALSE),

                            conditionalPanel("input.seqs_advanced",

                                numericInput("mv.count",
                                    "Minimum variant count",
                                    value=2, min=1, max=20, step=1),

#Explanation: The minimum variant count

                                p("Requires each mutation to occur in
                                the alignment at least this many times
                                for inclusion, e.g. 2 omits
                                singleton mutations.")
# from swarm selection

#                                 numericInput("mv.count",
#                                     "Maximum tolerated insertion length",
#                                     value=2, min=1, max=20, step=1),

#                                 numericInput("mv.count",
#                                     "Maximum tolerated deletion length",
#                                     value=2, min=1, max=20, step=1),

#                                p("To do: specify indel lengths that would be tolerated for inclusion."),
#                                p("Could also add support to include and exclude sequences by name.")

                            ), # end of advanced options conditional panel
hr(),
                            p(strong("Download Results")),
                            p(downloadButton('downloadSelectedSites',
                                "Table of selected sites")),

                            p(downloadButton('downloadSwarmsetNames',
                                "Selected sequence names")),

                            p(downloadButton('downloadSwarmsetSequences',
                                "Selected sequences in Raw/Table format")),

                            p(downloadButton('downloadFastaSwarmsetSequences',
                                "Selected sequences in Fasta format")),

                            p(downloadButton('downloadSwarmsetConcatamers',
                                "Selected concatamers in Raw/Table format")),

                            p(downloadButton('downloadMessages',
                                'Output messages')),

                            p(downloadButton('downloadSwarmsetLogos',
                                "Logos plot of selected concatamers")),

                            p("Options for Logos"),
                                radioButtons("logo.format", label=NA,
                                    choices = logo.formats, inline=T),
                                checkboxInput("logo.dotify",
                                    "Use blank for TF states", T),
                                checkboxInput("logo.stratify",
                                    "Make one plot per timepoint", T),
                                checkboxInput("logo.sort",
                                    "Sort sites from N- to C-term", T)
                        ) # end of select_seqs conditional panel
                    ),

                    mainPanel(
                        conditionalPanel("input.select_seqs",
                            code("Selected Sequences",
                                tableOutput("selectSequences")))
                    )
                )
            ), # end of tab panel

            tabPanel("About LASSIE", id="basics",
                includeMarkdown("basics.Rmd")
            ),  # end of tab panel
            tabPanel("Detailed Help", id="details",
                includeMarkdown("details.Rmd")
            ),  # end of tab panel
            tabPanel("License", id="license",
                includeMarkdown("legalese.Rmd")
            )  # end of tab panel
        ) # end of tabset
    ) # end of page
)

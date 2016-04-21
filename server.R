library(shiny)
#library(DT)

options(shiny.maxRequestSize=50 *1024^2) # 50MB
options(shiny.trace=F)

shinyServer(function(input, output, session) {

    S <- lassie::swarmtools()
    SS <- NULL
    infile <- NULL
#    values <- reactiveValues(infile = NULL)
#    values <- reactiveValues(cutoff = 0)

    updateInfile <- reactive({

#        gotInput()

        if (input$demo_data)
            infile = system.file("extdata", "CH505-gp160.fasta",
                package="lassie")

        if (!is.null(input$aas_file$datapath))
            infile = input$aas_file$datapath

        infile
    })

    gotInput <- reactive({

        validate(
            need(input$demo_data | input$aas_file$datapath | input$aas_file$name,
                    "Please select a protein alignment file."))
    })

    updateInputs <- reactive({

        validate(
	    need(input$tf_loss_when_up <= input$tf_loss_cutoff,
	    message = "Please review your settings.\n'Order sites by when they first reach this % TF loss' cannot exceed 'Exclude sites with % peak TF loss below'.",
        ))

        lassie::set.tf.loss.cutoff(
            lassie::set.alignment.file(updateParser(), updateInfile(),
                alignment_format=input$aln.format),
            as.numeric(input$tf_loss_cutoff), as.numeric(input$tf_loss_when_up))
    })

    updateSwarmset <- reactive({
        suppressMessages(lassie::swarmset(updateInputs(), 
	    min_counts=input$mv.count))
    })

    testParsing <- reactive({
        validate(
            need(
                all(grepl(input$tp.sep, rownames(S$aas_aln))),
                message = paste0("Please choose a different field separator; ",
                    input$tp.sep, " does not occur in all sequence names.")),
            need(
                all(sapply(1:nrow(S$aas_aln), function(i)
                    length(unlist(strsplit(rownames(S$aas_aln)[i],
                        input$tp.sep, fixed=T)))) >= input$tp.pos),
                message = "Please review your sequence names or choose a lower value for sample timepoint field.  Not all sequence names have this many fields.")
        )
    })

    updateParser <- reactive({

        S$pngs2o=input$map.pngs2o

        S <- lassie::set.tp.parser(S,
            lassie::create.timepoint.parser(field.sep =
                as.character(input$tp.sep),
                field.num = as.numeric(input$tp.pos)))

#        testParsing()

        return ( S )
    })

    output$nSelectedSites <- renderText({

#        S <- updateInputs()

        if (!is.null(S$aas_aln) & !is.null(S$selected_sites))
            paste0("The table below summarizes sites with no less than ",
                input$tf_loss_cutoff, "% TF loss.  There are ",
                nrow(S$selected_sites), " such sites.")# among ",
#                ncol(S$aas_aln), " aligned columns.")
    })

    observe({
        tf_cutoff <- input$tf_loss_cutoff
    	updateSliderInput(session, "tf_loss_when_up", max=tf_cutoff)
    })

    output$tabulateByTimepoint <- renderTable({
	validate(
            need(input$demo_data | !is.null(input$aas_file$datapath),
            "")
        )

	S <- updateInputs()

# S$n_per_timepoint is a named one-D vector i.e. table object
# need to convert it to a data frame.  
# transpose function is used to list all in one named row.
	try ( as.data.frame(t(as.data.frame.vector(S$n_per_timepoint, optional=T)), row.names=c('n')) )
    })

#    output$selectedSites <- DT::renderDataTable({
    output$selectedSites <- renderDataTable({

	validate(
            need(input$demo_data | !is.null(input$aas_file$datapath),
            "Please select an alignment file")
        )

        S <- updateInputs()

        if (!is.null(S$selected_sites)) {
            if (nrow(S$selected_sites) > 0) {
                if (!all(grepl(input$tp.sep, rownames(S$aas_aln))))
                    as.data.frame(paste0("Please choose a different field separator; ",
                        input$tp.sep, " does not occur in all sequence names."))
                    else if (!all(sapply(1:nrow(S$aas_aln), function(i)
                        length(unlist(strsplit(rownames(S$aas_aln)[i],
                            input$tp.sep, fixed=T)))) >= input$tp.pos))
                        as.data.frame("Please review your sequence names or choose a lower value for sample timepoint field.  Not all sequence names have this many fields.")
                return ( S$selected_sites )
            }
        }
    })

    selectSequences <- reactive({
        validate(
            need(input$demo_data | !is.null(input$aas_file$datapath), # | input$aas_file$name,
                "Please check the 'Send selected sites to Step 2' box in Step 1."))
    })

    output$selectSequences <- renderTable({

#        validate(need(input$demo_data | input$aas_file$name != "",
#            "Please check the 'Send selected sites to Step 2' box in Step 1."))
#        selectSequences()
        validate(
            need(input$demo_data | !is.null(input$aas_file$datapath), # | input$aas_file$name,
                "Please check the 'Send selected sites to Step 2' box in Step 1."))
#        if (input$select_seqs &
#                (input$demo_data | !is.null(input$aas_file$datapath)))
            lassie::tabulate.swarmset(updateSwarmset())
#        else
#            NULL
    })

    output$showMessages <- renderTable({
#        if (input$select_seqs &
#                (input$demo_data | !is.null(input$aas_file$datapath))) {
selectSequences()
            SS <- updateSwarmset()
            message = SS$working_swarm$message
            as.data.frame(message, row.names=NULL)
#        } else {
#            return (NULL)
#        }
    })

    output$showSummary <- renderText({

#        if (input$select_seqs &
#                (input$demo_data | !is.null(input$aas_file$datapath))) {
selectSequences()
            SS <- updateSwarmset()

            paste("Selected", nrow(lassie::tabulate.swarmset(SS)),
                #"of", nrow(S$aas_aln),
                "sequences to represent selected sites.")
#        } else {
#            c("Please select sites in the 'Step 1' tab.")
#        }
     })

    output$downloadSelectedSites <- downloadHandler(
        filename = function() {
            paste0('lassie-sitelist-', Sys.Date(), '.txt')
        },
        content = function(file) {
            write.table(print(updateInputs()), file, quote=F, row.names=F)
        }
    )

    output$downloadSwarmsetNames <- downloadHandler(
        filename = function() {
            paste0('lassie-seqnames-', Sys.Date(), '.txt')
        },
        content = function(file) {
            write.table(capture.output(summary(updateSwarmset())),
                file, quote=F, col.names=F, row.names=F)
        }
    )

    output$downloadSwarmsetSequences <- downloadHandler(
        filename = function() {
            paste0('lassie-sequences-', Sys.Date(), '.txt')
        },
        content = function(file) {
            write.table(lassie::export.swarmset(updateSwarmset()),
                file, quote=F, col.names=F, row.names=T)
        }
    )

    output$downloadFastaSwarmsetSequences <- downloadHandler(
        filename = function() {
            paste0('lassie-sequences-', Sys.Date(), '.fasta')
        },
        content = function(file) {
            lassie::export.fasta.swarmset(updateSwarmset(), file)
        }
    )

    output$downloadSwarmsetConcatamers <- downloadHandler(
        filename = function() {
            paste0('lassie-concatamers-', Sys.Date(), '.txt')
        },
        content = function(file) {
            write.table(lassie::tabulate.swarmset(updateSwarmset()),
                file, quote=F, col.names=F)
        }
    )

    output$downloadSwarmsetLogos <- downloadHandler(

        filename = function() {
            paste0('lassie-swarmset-logos-', Sys.Date(), '.',
                ifelse(input$logo.stratify, "zip", input$logo.format))
        },
        content = function(file) {
            tmp.file = lassie::make.timepoint.logos(updateSwarmset(),
                sort_stacks=input$logo.sort,
                stratify=input$logo.stratify,
                dotify=input$logo.dotify,
                logo_format = input$logo.format)

            img.filenames <- basename(tmp.file)
            img.dir <- dirname(tmp.file)

            if (length(tmp.file) > 1 | input$logo.stratify) {
                curr.dir <- getwd()
                setwd(img.dir)
                zip(file, img.filenames)
                setwd(curr.dir)
            } else {
                file.rename(tmp.file, file)
            }
        }
    )

    output$downloadMessages <- downloadHandler(
        filename = function() {
            paste0('lassie-messages-', Sys.Date(), '.txt')
        },
        content = function(file) {
#            if (input$select_seqs &
#                (input$demo_data | !is.null(input$aas_file$datapath))) {
            selectSequences()

                SS <- updateSwarmset()
                message = SS$working_swarm$message

                write.table(as.data.frame(message, row.names=NULL),
                    file, quote=F, col.names=F, row.names=F)
#            } else {
#                return (NULL)
#            }
        }
    )

})

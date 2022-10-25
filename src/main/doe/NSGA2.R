#title: NSGA2
#help: Minimizes a multidimensional function to approximate its Pareto front and Pareto set
#type: Optimization
#authors: Yann Richet <yann.richet@irsn.fr>; Heike Trautmann <trautmann@statistik.uni-dortmund.de>, Detlef Steuer <steuer@hsu-hamburg.de> and Olaf Mersmann <olafm@statistik.uni-dortmund.de>
#references: Deb, K., Pratap, A., and Agarwal, S.. A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II. IEEE Transactions on Evolutionary Computation, 6 (8) (2002), 182-197.
#require: future; templr; mco
#options: output_index='1';constraints_index='';popsize_quarter='10';generations = '10';cprob = '0.7';cdist = '5';mprob = '0.2';mdist = '10'; seed='123'
#options.help: output_index=Main output indexes;constraints_index=Output constraints indexes;popsize_quarter=Size of population div. by 4;generations=Number of generations to breed;cprob=Crossover probability;cdist=Crossover distribution index;mprob=Mutation probability;mdist=Mutation distribution index

NSGA2 <- function(options) {
    library(mco)
    
    algorithm = new.env()

    algorithm$output_index <- as.integer(options$output_index)
    if (is.character(options$constraints_index)) { # remove trailing quotes
        options$constraints_index=gsub("'","",options$constraints_index,fixed = T)
        options$constraints_index=gsub('"',"",options$constraints_index,fixed = T)
    }
    algorithm$constraints_index <- as.integer(options$constraints_index)

    algorithm$popsize <- 4*as.integer(options$popsize_quarter)
    algorithm$generations <- as.integer(options$generations)
    algorithm$cprob <- as.numeric(options$cprob)
    algorithm$cdist <- as.numeric(options$cdist)
    algorithm$mprob <- as.numeric(options$mprob)
    algorithm$mdist <- as.integer(options$mdist)
        
    algorithm$seed <- as.integer(options$seed)
    
    algorithm$id <- floor(1000*runif(1))
    
    return(algorithm)
}

getInitialDesign <- function(algorithm, input, output) {
    if (length(algorithm$output_index)+length(algorithm$constraints_index) != length(output)) {
        algorithm$output_index = 1:length(output)
        algorithm$output_index[algorithm$constraints_index] <- NA
        algorithm$output_index = na.omit(algorithm$output_index)

        algorithm$constraints_index = (1:length(output))[-algorithm$output_index]
    }

    set.seed(algorithm$seed)
    algorithm$input <- input
    algorithm$output <- output  
    algorithm$d <- length(input)
    library(future)
    library(templr)
    wd = getwd()
    algorithm$job = future(evaluator=plan("multisession"),lazy = FALSE,{
        sink(file.path(wd,paste0('NSGA2_',algorithm$id,'.out')),type='output')
        print("Starting nsga2()")
        set.seed(algorithm$seed)
        if (length(algorithm$constraints_index)>0) {
            o = nsga2(fn = function(x) {ask_Y(id=algorithm$id,matrix(x,ncol=algorithm$d))[,algorithm$output_index]},
                      constraints = function(x) {ask_Y(id=algorithm$id,matrix(x,ncol=algorithm$d))[,algorithm$constraints_index]},
                      idim = algorithm$d, odim=length(algorithm$output_index), cdim=length(algorithm$constraints_index),
                      popsize =algorithm$popsize, generations = algorithm$generations, 
                      cprob = algorithm$cprob,cdist = algorithm$cdist,mprob = algorithm$mprob,mdist = algorithm$mdist,vectorized = TRUE, 
                      lower.bounds=min_input(input), upper.bounds=max_input(input))
        } else {
            o = nsga2(fn = function(x) {ask_Y(id=algorithm$id,matrix(x,ncol=algorithm$d))},
                      idim = algorithm$d, odim=length(algorithm$output_index),
                      popsize =algorithm$popsize, generations = algorithm$generations, 
                      cprob = algorithm$cprob,cdist = algorithm$cdist,mprob = algorithm$mprob,mdist = algorithm$mdist,vectorized = TRUE, 
                      lower.bounds=min_input(algorithm$input), upper.bounds=max_input(algorithm$input))
        }
        print("nsga2() ended")
        if (is.null(o)) 
          print(traceback())
        else {
            saveRDS(o,file=paste0("NSGA2_",algorithm$id,".rds")) # serialize optim result for later displayResults...
        }
        sink(type='output')
        ask_Y(id=algorithm$id,matrix(NaN,ncol=algorithm$d))
    })
    algorithm$i = 0
    
    Sys.sleep(.1)
    
    Xn = ask_X(id=algorithm$id)
    algorithm$s = nrow(Xn)

    colnames(Xn) <- names(algorithm$input)
    return(Xn)
}

getNextDesign <- function(algorithm, X, Y) {
    algorithm$i = algorithm$i + 1
    
    y = Y[(nrow(Y)-algorithm$s+1):nrow(Y),]
    if (all(is.na(y))) {
        tell_Y(id=algorithm$id,NULL)
        return(NULL)
    } else tell_Y(id=algorithm$id,y)
    
    Sys.sleep(.1)

    Xn = ask_X(id=algorithm$id)
    if (is.null(Xn)) return(NULL)
    algorithm$s = nrow(Xn)

    colnames(Xn) <- names(algorithm$input)
    return(Xn)
}

displayResults <- function(algorithm, X, Y) {
    algorithm$files <- paste0("plot[",algorithm$i,"].png")
    png(file = algorithm$files, height = 600, width = 600)
    Yo = apply(Y[,algorithm$output_index,drop=F],1,prod)
    Yc = apply(Y[,algorithm$constraints_index,drop=F],1,function(yc) all(yc>0) )
    red=(Yo-min(Yo))/(max(Yo)-min(Yo))
    pairs(cbind(X,Y),col=rgb(r=red,g=0,b=1-red, alpha = 0.1+0.9*Yc))
    dev.off()
    
    html = paste(sep = "",
                 "<HTML name='points'>",
                 "<img src='",algorithm$files,"' width='600' height='600'/>",
                 "<br/></HTML>")

    if (file.exists(paste0("NSGA2_",algorithm$id,".rds"))) {
        algorithm$optim = readRDS(paste0("NSGA2_",algorithm$id,".rds"))
        sol = which(algorithm$optim$pareto.optimal)
        if (length(sol)>0) html=paste0(html,
           "<min>[",paste0(collapse = ",",algorithm$optim$value[sol]),"]</min>",
           "<argmin>",paste0(collapse = ",\n","[",paste0(collapse = ",",algorithm$optim$par[sol,]),"]"),"</argmin>")
    }

    ## return HTML string containing plot image
    return(html)
}

displayResultsTmp <- displayResults

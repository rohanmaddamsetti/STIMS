"""
STIMS.jl by Rohan Maddamsetti.

Recommended usage is to run the RunSTIMS function from the Julia REPL.

For examples for how to run from the Julia REPL, see the function run_examples().

Examples for how to run from the command-line:

julia STIMS.jl ../results/LTEE-metagenome-mutations.csv ../results/REL606_IDs.csv ../data/neutral_compilation.csv -o ../results/gene-modules/STIMS-jl-test-figures/STIMS-plot.pdf

julia STIMS.jl ../results/SLiM-5000gen-v03.csv ../results/SLiM_geneIDs.csv ../results/SLiM_test_gene_module.csv -o ../results/gene-modules/STIMS-jl-test-figures/STIMS-plot.pdf

ISSUES:

1) crashes when the data is too huge, due to plotting too many points.
   Nkrumah proposed the following solution:
   When the data reaches a critical size, plot a random sample of the
   data to reduce it to some size.

"""

module STIMS

## exported interface-- the only functions that users should use.
export RunSTIMS, RunSTIMS_on_data

using DataFrames, DataFramesMeta, CSV, StatsBase, FLoops, RCall, ArgParse

## See this blog post about using ggplot within Julia:
## https://avt.im/blog/2018/03/23/R-packages-ggplot-in-julia
@rlibrary ggplot2
@rlibrary scales


function make_cumulative_muts_pop_df(gene_module_mutations_df, pop, final_time)
    gene_module_mutations_in_pop_df = @rsubset(gene_module_mutations_df,
                                               :Population == pop)
    
    if (nrow(gene_module_mutations_in_pop_df) == 0)
        ## no mutations in this pop
        almost_done_df = DataFrame(Population = pop,
                                   t0 = final_time,
                                   count = 0,
                                   cs = 0)
    else
        summary_df = @chain gene_module_mutations_in_pop_df begin
            groupby([:Population, :t0])
            combine(nrow => :count)
            sort(:t0)
            transform(:count => cumsum => :cs)
        end
        
        ## since the final_time is not in ret_df, add one final row,
        ## for nicer plots.
        final_row_df = DataFrame(Population = pop,
                                 t0 = final_time,
                                 count = 0,
                                 cs = maximum(summary_df.cs))
        
        almost_done_df = vcat(summary_df, final_row_df)
    end
    ## add an row for Time == 0 (for nicer plots).
    init_row_df = DataFrame(Population = pop,
                            t0 = 0,
                            count = 0,
                            cs = 0)
    
    pop_df = vcat(init_row_df, almost_done_df)
    return pop_df
end


function calc_cumulative_muts(gene_module_df,
                              gene_mutation_data, genome_metadata, pop_level_vec)
    ## look at accumulation of stars over time
    ## in other words, look at the rates at which the mutations occur over time.
    ## To normalize, we need to supply the number of sites at risk
    ## (such as sum of gene length).

    gene_module_mut_data = @rsubset(gene_mutation_data, :Gene in gene_module_df.Gene)
    gene_module_metadata = @rsubset(genome_metadata, :Gene in gene_module_df.Gene)
    
    ## normalize by the total length of genes
    ## in the given module (in d.metadata).
    my_genes = @chain gene_module_metadata begin
        select(:Gene, :gene_length)
        unique([:Gene, :gene_length])
    end
    
    normalization_constant = sum(my_genes.gene_length)
    
    final_time = maximum(gene_mutation_data.t0) + 1
    ## + 1 so that final_time is outside of the actual data.
    
    ## we are going to concatenate the summarized data to this empty DataFrame.
    c_dat = DataFrame()
    
    for pop in pop_level_vec
        pop_df = make_cumulative_muts_pop_df(gene_module_mut_data, pop, final_time)
        ## concatenate the pop_df to cumulative_df.
        c_dat = vcat(c_dat, pop_df)
    end
    
    transform!(c_dat, :cs => ByRow(cs -> cs/normalization_constant) => :normalized_cs)
    ## possible TODO: remove any NA values that arise?
    
    return c_dat
end


function generate_cumulative_mut_subset(gene_mutation_data, genome_metadata,
                                        pop_level_vec,
                                        subset_size)
    ## This function calculates cumulative mutations for a random gene set.
    rando_genes = sample(genome_metadata.Gene, subset_size; replace=false)
    rando_genes_df = @rsubset(genome_metadata, :Gene in rando_genes)
    
    c_mut_subset = calc_cumulative_muts(rando_genes_df, gene_mutation_data, genome_metadata, pop_level_vec)
    return c_mut_subset
end


function calc_traj_pvals(gene_module_df,
                         gene_mutation_data, genome_metadata,
                         pop_level_vec; N = 10000)    
    #= calculate the tail probabilities of the true cumulative mutation trajectory
    of a given vector of genes (a 'module'), based on resampling
    random sets of genes. Returns the upper tail of null distribution,
    or P(random trajectory >= the actual trajectory).
    Output: a dataframe with three columns: Population, count, p.val.
    =#

    nthreads = Threads.nthreads() ## get the number of threads available to Julia.
    
    ## each sample has the same cardinality as gene_module_df.Gene.
    subset_size = length(gene_module_df.Gene)

    data_trajectory = calc_cumulative_muts(gene_module_df,
                                           gene_mutation_data,
                                           genome_metadata,
                                           pop_level_vec)
    
    data_final_normalized_cs_summary = @chain data_trajectory begin
        groupby(:Population)
        combine(:normalized_cs => maximum => :data_final_normalized_cs)
    end

    @floop ThreadedEx(basesize = N ÷ nthreads) for _ in 1:N 
        randomized_trajectory = generate_cumulative_mut_subset(gene_mutation_data,
                                                               genome_metadata,
                                                               pop_level_vec,
                                                               subset_size)
        
        randomized_final_normalized_cs_summary = @chain randomized_trajectory begin
            groupby(:Population)
            combine(:normalized_cs => maximum => :randomized_final_normalized_cs)
            ## compare the randomized subset to the actual data:
            innerjoin(data_final_normalized_cs_summary, on = :Population)
            combine(:, [:randomized_final_normalized_cs, :data_final_normalized_cs] =>
                    ByRow((a, b) -> a > b) => :greater_than_data)
        end
        
        greater_than_data_vec = randomized_final_normalized_cs_summary.greater_than_data
        @reduce() do (count_vec = zeros(length(pop_level_vec)); greater_than_data_vec)
            ## now update the counts.
            count_vec .+= greater_than_data_vec
        end
    end
    
    uppertail_prob_df = DataFrame("Population" => pop_level_vec,
                                  "count" => count_vec,
                                  "pval" => count_vec/N); 
    return uppertail_prob_df
end


function get_middle_trajectories(bootstrapped_trajectories, alphaval, N)
    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories
    ## from each population, for a two-sided test.
    ## usual default is alphaval == 0.05.
    
    trajectory_summary = @chain bootstrapped_trajectories begin
        groupby([:bootstrap_replicate, :Population])
        combine(:normalized_cs => maximum => :final_norm_cs)
    end

    num_traj_to_slice = Integer(N * alphaval) ÷ 2
    
    top_trajectories = @chain trajectory_summary begin
        sort([:Population, :final_norm_cs])
        groupby(:Population)
        combine(x -> last(x, num_traj_to_slice))
        select(:Population, :bootstrap_replicate)
        @rtransform(:in_top = true)
    end

    bottom_trajectories = @chain trajectory_summary begin
        sort([:Population, :final_norm_cs])
        groupby(:Population)
        combine(x -> first(x, num_traj_to_slice))
        select(:Population, :bootstrap_replicate)
        @rtransform(:in_bottom = true)
    end

    middle_trajectories = @chain bootstrapped_trajectories begin
        leftjoin(top_trajectories, on = [:Population, :bootstrap_replicate])
        leftjoin(bottom_trajectories, on = [:Population, :bootstrap_replicate])
        @rsubset(ismissing(:in_top))
        @rsubset(ismissing(:in_bottom))
        select(Not(:in_top))
        select(Not(:in_bottom))
    end
    
    return middle_trajectories
end


function plot_base_layer(gene_mutation_data, genome_metadata, pop_level_vec;    
                         subset_size = 50, N = 1000, alphaval = 0.05,
                         my_color = "gray")
    #= This plot visualizes a two-tailed test (alphaval = 0.05)
    against a bootstrapped null distribution.
    Throughout, plots use the minimum subsample size to subsample
    the null distribution, to increase the variance in order to
    make a conservative comparison. 
    =#

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    bootstrapped_trajectories = DataFrame()
    for bootstrap_rep in 1:N 
        randomized_traj = generate_cumulative_mut_subset(gene_mutation_data,
                                                         genome_metadata,
                                                         pop_level_vec,
                                                         subset_size)
        ## add a column for bootstrap replicate
        @transform!(randomized_traj, :bootstrap_replicate = bootstrap_rep)
        bootstrapped_trajectories = vcat(bootstrapped_trajectories, randomized_traj)
    end

    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories
    ## from each population, for a two-sided test.
    ## default is alphaval == 0.05.
    middle_trajs = get_middle_trajectories(bootstrapped_trajectories, alphaval, N)
    ## for better plotting, divide x-axis labels by 1000.
    transform!(middle_trajs, :t0 => ByRow(t -> t/1000) => :Time)

    ## R function for plotting better y-axis labels.
    ## see solution here for nice scientific notation on axes.
    ## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales

    fancy_scientific = R"""function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))}"""

    p = ggplot(middle_trajs, aes(x=:Time, y=:normalized_cs)) +
        ylab("Cumulative mutations (normalized)") +
        theme_classic() +
        geom_point(size=0.2, color=my_color) +
        theme(var"axis.title.x" = element_text(size=13),
              var"axis.title.y" = element_text(size=13),
              var"axis.text.x"  = element_text(size=13),
              var"axis.text.y"  = element_text(size=13)) +
                  scale_y_continuous(labels=fancy_scientific,
                                     breaks = R"scales::extended_breaks(n = 6)",
                                     limits = R"c(0, NA)") +
                  facet_wrap(R".~Population", scales="free", nrow=4) +
                  xlab("Generations (x 1,000)")
    return p
end


function add_cumulative_mut_layer(p, layer_df; my_color="black")
    ## take a ggplot object output by plot_cumulative_muts, and add an extra layer.
    p = p +
        geom_point(data = layer_df,
                   aes(x = :Time, y = :normalized_cs),
                   color = my_color, size = 0.2) +
                       geom_step(data = layer_df, aes(x = :Time, y = :normalized_cs),
                                 size = 0.2, color = my_color)
    return p
end


function run_LTEE_analyses()
    ## This function is for checking correctness,
    ## by comparison to my STIMS implementation in R.
    
    ## get the lengths of all protein-coding genes in REL606.
    ## This excludes genes in repetitive regions of the genome.
    ## See Section 4.3.1 "Removing mutations in repetitive regions of the genome"
    ## in Ben Good's LTEE metagenomics paper for more details.
    ## This filtering is done in my python script printEcoliIDs.py.
    ##Do by running:
    ## python printEcoliIDs.py -i ../data/REL606.7.gbk --mask ../data/REL606.L20.G15.P0.M35.RM-edited.mask.gd > ../results/REL606_IDs.csv
    REL606_genes = CSV.read("../results/REL606_IDs.csv", DataFrame) 
    
    mutation_data = CSV.read("../results/LTEE-metagenome-mutations.csv", DataFrame)
    
    ## for documentation of this syntax, see:
    ## https://juliadata.github.io/DataFramesMeta.jl/stable/#Comparison-with-dplyr-and-LINQ
    
    gene_mutation_data = @chain mutation_data begin
        ## important: innerjoin only includes genes that pass REL606 annotation
        ## filters. In particular, genes that overlap with masked regions are
        ## excluded, even if the metagenomic data included mutations that
        ## are called in non-repetitive regions of those genes.
        transform(:t0 => ByRow(t -> t/1000) => :Time)
        innerjoin(REL606_genes, on = :Gene)
        @rsubset(:Gene != "intergenic")
    end
    
    
    ## This is only useful for LTEE data.
    ## This constant is to make sure that all pops are in the levels
    ## of the Population factor after mergers, etc.
    LTEE_pop_level_vec =  ["Ara-5","Ara-6", "Ara+1", "Ara+2",
                           "Ara+4", "Ara+5", "Ara-1", "Ara-2",
                           "Ara-3", "Ara-4", "Ara+3", "Ara+6"]
    
    
    ## get mutation parallelism in the LTEE genomes published in Tenaillon (2016).
    ## This is used for filtering essential genes
    ## in the purifying selection control analysis,
    ## and is used in the positive selection control analysis.
    nonmut_genomics = CSV.read("../data/tenaillon2016-nonmutator-parallelism.csv",
                               DataFrame, normalizenames=true)
    ## rename the :Gene_name column to :Gene.
    rename!(nonmut_genomics,:Gene_name => :Gene)
    ## make sure these genes passed the filters on REL606.genes.
    @rsubset!(nonmut_genomics, :Gene in REL606_genes.Gene)
    sort!(nonmut_genomics, :G_score, rev = true)
    
    hypermut_genomics = CSV.read("../data/tenaillon2016-mutator-parallelism.csv",
                                 DataFrame, normalizenames=true)
        ## rename the :Gene_name column to :Gene.
    rename!(hypermut_genomics,:Gene_name => :Gene)
    ## make sure these genes passed the filters on REL606.genes.
    @rsubset(hypermut_genomics, :Gene in REL606_genes.Gene)    
    sort!(hypermut_genomics, :G_score, rev = true)
    
    top_nonmut_genomics = first(nonmut_genomics, 50)
    top_hypermut_genomics = first(hypermut_genomics, 50)
    
    neutral_genes = CSV.read("../data/neutral_compilation.csv", DataFrame)
    ## make sure that only loci that pass filters are included in the analysis.
    @rsubset!(neutral_genes, :Gene in REL606_genes.Gene)
    
    c_neutral_genes = calc_cumulative_muts(neutral_genes,
                                           gene_mutation_data,
                                           REL606_genes,
                                           LTEE_pop_level_vec)
    ## for better plotting, divide x-axis labels by 1000.
    transform!(c_neutral_genes, :t0 => ByRow(t -> t/1000) => :Time)


    neutral_pvals = calc_traj_pvals(neutral_genes,
                                    gene_mutation_data, REL606_genes,
                                    LTEE_pop_level_vec)
    println("gold standard relaxed selection results:")
    println(neutral_pvals) ## print the results

    
    neutral_base_layer = plot_base_layer(
        gene_mutation_data, REL606_genes, LTEE_pop_level_vec,
        subset_size = length(unique(neutral_genes.Gene)))
    
    ## Figure 4: plot of "gold standard" neutral genes.
    Fig4 = add_cumulative_mut_layer(
        neutral_base_layer,
        c_neutral_genes)
    ggsave("../results/gene-modules/STIMS-jl-test-figures/Fig4.pdf", Fig4)

    
    essential_genes = CSV.read("../data/Couce2017-LTEE-essential.csv", DataFrame)
    essential_genes = @chain essential_genes begin
        innerjoin(REL606_genes, on = :Gene)
        @rsubset(!ismissing(:locus_tag))
    end
    ## a significant proportion of genes under positive selection in the LTEE are
    ## essential genes, as reported in Maddamsetti et al. (2017).
    ## filter these ones out, since we are interested in purifying selection.
    ## 21 out of 50 top non-mut genes are essential.
    nonmut_top_hit_essential = @rsubset(essential_genes,
                                        :Gene in top_nonmut_genomics.Gene)
    ## what about the hypermutators? 3 out of 50 top hypermut genes.
    hypermut_top_hit_essential = @rsubset(essential_genes, :Gene in top_hypermut_genomics.Gene)
    
    ## filtering out top G-score genes in the LTEE genomics dataset.
    purifying_genes = @chain essential_genes begin
        @rsubset(!(:Gene in nonmut_top_hit_essential.Gene))
        @rsubset(!(:Gene in hypermut_top_hit_essential.Gene))
    end

c_purifying_genes = calc_cumulative_muts(purifying_genes,
                                         gene_mutation_data,
                                         REL606_genes, LTEE_pop_level_vec)
## for better plotting, divide x-axis labels by 1000.
transform!(c_purifying_genes, :t0 => ByRow(t -> t/1000) => :Time)

purifying_base_layer = plot_base_layer(
    gene_mutation_data, REL606_genes, LTEE_pop_level_vec,
    subset_size=length(unique(purifying_genes.Gene)))

##  Yes. evidence of purifying selection in these genes based on my test.
Fig5 = add_cumulative_mut_layer(purifying_base_layer,
                                c_purifying_genes)
ggsave("../results/gene-modules/STIMS-jl-test-figures/Fig5.pdf", Fig5)


## calculate more rigorous statistics than the figures.
purifying_pvals = calc_traj_pvals(purifying_genes,
                                  gene_mutation_data, REL606_genes,
                                  LTEE_pop_level_vec)
println("gold standard purifying selection results:")
println(purifying_pvals) ## print the results


## now look at positive selection genes.
rando_plot = plot_base_layer(gene_mutation_data, REL606_genes,
                             LTEE_pop_level_vec)

## 1) plot top genes in non-mutators.
c_top_nonmuts =  calc_cumulative_muts(top_nonmut_genomics,
                                      gene_mutation_data,
                                      REL606_genes, LTEE_pop_level_vec)
## for better plotting, divide x-axis labels by 1000.
transform!(c_top_nonmuts, :t0 => ByRow(t -> t/1000) => :Time)


Fig6 = add_cumulative_mut_layer(rando_plot, c_top_nonmuts)
ggsave("../results/gene-modules/STIMS-jl-test-figures/Fig6.pdf",Fig6)

## calculate more rigorous statistics than the figures.
top_nonmut_pvals = calc_traj_pvals(top_nonmut_genomics,
                                   gene_mutation_data,
                                   REL606_genes,
                                   LTEE_pop_level_vec)
println("gold standard positive selection results:")
println(top_nonmut_pvals) ## print the results
end


function RunSTIMS(mutation_csv_path, genome_metadata_csv_path,
                  genelist_csv_path, outfile = "STIMS-plot.pdf")

    mutation_data = CSV.read(mutation_csv_path, DataFrame)
    genome_metadata = CSV.read(genome_metadata_csv_path, DataFrame) 
   
    ## for documentation of this syntax, see:
    ## https://juliadata.github.io/DataFramesMeta.jl/stable/#Comparison-with-dplyr-and-LINQ
    
    gene_mutation_data = @chain mutation_data begin
        transform(:t0 => ByRow(t -> t/1000) => :Time)
        innerjoin(genome_metadata, on = :Gene)
        @rsubset(:Gene != "intergenic")
    end
    
    pop_level_vec = unique(gene_mutation_data.Population)
    
    gene_module_df = CSV.read(genelist_csv_path, DataFrame)
    ## make sure that only loci in genome_metadata are analyzed.
    @rsubset!(gene_module_df, :Gene in genome_metadata.Gene)

    c_gene_module = calc_cumulative_muts(gene_module_df,
                                         gene_mutation_data,
                                         genome_metadata,
                                         pop_level_vec)
    ## for better plotting, divide x-axis labels by 1000.
    transform!(c_gene_module, :t0 => ByRow(t -> t/1000) => :Time)
    
    pvals = calc_traj_pvals(gene_module_df,
                            gene_mutation_data,
                            genome_metadata,
                            pop_level_vec)

    println(pvals) ## print the results
    
    base_layer = plot_base_layer(
        gene_mutation_data, genome_metadata, pop_level_vec,
        subset_size = length(unique(gene_module_df.Gene)))
    
    Fig = add_cumulative_mut_layer(
        base_layer,
        c_gene_module)
    ## default figure dimensions.
    fig_height = 7
    fig_width = 7
    ## if plotting a single population, then plot a smaller figure.
    if (length(pop_level_vec) == 1)
        fig_height = 3
        fig_width = 3
    end
    ggsave(outfile, Fig, height = fig_height, width = fig_width)
end


## take dataframes directly as input, and don't make the plot.
function RunSTIMS_on_data(mutation_data, genome_metadata, gene_module_df)
   
    gene_mutation_data = @chain mutation_data begin
        transform(:t0 => ByRow(t -> t/1000) => :Time)
        innerjoin(genome_metadata, on = :Gene)
        @rsubset(:Gene != "intergenic")
    end
    
    pop_level_vec = unique(gene_mutation_data.Population)
    
    ## make sure that only loci in genome_metadata are analyzed.
    @rsubset!(gene_module_df, :Gene in genome_metadata.Gene)

    c_gene_module = calc_cumulative_muts(gene_module_df,
                                         gene_mutation_data,
                                         genome_metadata,
                                         pop_level_vec)

    pvals = calc_traj_pvals(gene_module_df,
                            gene_mutation_data,
                            genome_metadata,
                            pop_level_vec)
    return(pvals)
end


function run_examples()
    RunSTIMS("../results/LTEE-metagenome-mutations.csv", "../results/REL606_IDs.csv", "../data/neutral_compilation.csv", "../results/gene-modules/STIMS-jl-test-figures/neutral.pdf")
end


function run_SLiM_tests()
    RunSTIMS("../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator.csv",
             "../results/SLiM-results/SLiM_geneIDs.csv",
             "../results/SLiM-results/SLiM_neutral_module.csv",
             "../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator-neutral.pdf")

    RunSTIMS("../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator.csv",
             "../results/SLiM-results/SLiM_geneIDs.csv",
             "../results/SLiM-results/SLiM_positive_module.csv",
             "../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator-positive.pdf")

    RunSTIMS("../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator.csv",
             "../results/SLiM-results/SLiM_geneIDs.csv",
             "../results/SLiM-results/SLiM_purifying_module.csv",
             "../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator-purifying.pdf")

    RunSTIMS("../results/SLiM-results/SLiM-5000gen-OnePercent-Nonmutator.csv",
             "../results/SLiM-results/SLiM_geneIDs.csv",
             "../results/SLiM-results/SLiM_neutral_module.csv",
             "../results/SLiM-results/SLiM-5000gen-OnePercent-Nonmutator-neutral.pdf")

    RunSTIMS("../results/SLiM-results/SLiM-5000gen-OnePercent-Nonmutator.csv",
             "../results/SLiM-results/SLiM_geneIDs.csv",
             "../results/SLiM-results/SLiM_positive_module.csv",
             "../results/SLiM-results/SLiM-5000gen-OnePercent-Nonmutator-positive.pdf")

    RunSTIMS("../results/SLiM-results/SLiM-5000gen-OnePercent-Nonmutator.csv",
             "../results/SLiM-results/SLiM_geneIDs.csv",
             "../results/SLiM-results/SLiM_purifying_module.csv",
             "../results/SLiM-results/SLiM-5000gen-OnePercent-Nonmutator-purifying.pdf")
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "mutation_csv_path"
        help = "csv file of mutations, formatted for STIMS."
        arg_type = String
        required = true
        "genome_metadata_csv_path"
        help = "csv file of genes in the genome, formatted for STIMS."
        arg_type = String
        required = true
        "gene_set_path"
        help = "a text file starting with the column name 'Gene', with one gene per line. This contains the set of genes that will be analyzed by STIMS."
        arg_type = String
        required = true
        "--outfile", "-o"
        help = "path for the figure. default is STIMS-plot.pdf."
        arg_type = String
        default = "STIMS-plot.pdf"
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    mutation_csv_path = parsed_args["mutation_csv_path"]
    genome_metadata_csv_path = parsed_args["genome_metadata_csv_path"]
    gene_set_csv_path = parsed_args["gene_set_path"]
    outfile = parsed_args["outfile"]
    
    RunSTIMS(mutation_csv_path, genome_metadata_csv_path,
             gene_set_csv_path, outfile)
end


## only run main() when STIMS.jl is invoked from the command-line.
if length(ARGS) > 0
    main()
end


end ## end of module STIMS.

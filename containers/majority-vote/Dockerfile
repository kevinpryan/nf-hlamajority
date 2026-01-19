FROM kevinr9525/r-basic:dev

COPY parse_outputs_majority_vote.R /bin/
COPY majority_vote_combine_nd_ad_lens.R /bin/ 

RUN chmod +x /bin/majority_vote_combine_nd_ad_lens.R && \
    chmod +x /bin/parse_outputs_majority_vote.R

RUN Rscript -e "install.packages(c('tidytable'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
COPY lib/* ./bin/lib/

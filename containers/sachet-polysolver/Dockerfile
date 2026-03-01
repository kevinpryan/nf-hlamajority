FROM sachet/polysolver:v4
COPY shell_call_hla_type /home/polysolver/scripts/shell_call_hla_type
COPY changes.sh .
RUN chmod +x /home/polysolver/scripts/shell_call_hla_type && chmod +x changes.sh 
CMD bash changes.sh

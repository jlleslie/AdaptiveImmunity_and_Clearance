# AdaptiveImmuinity_and_Clearance
Data and results interrogating the role of adaptive immunity in clearing C. difficile infection.

Overview
--------
    |- README          # Overview of all content
    |
    |- data            # Raw and primary data
    |  |- Adaptiveimmuneclear_metadata_noD40.42.txt         # Metadata for mouse experimental groups
    |  |- CDIclear.final.0.03.cons.taxonomy  # Raw taxonomy file from mothur analysis
    |  |- clearance.formatted.taxonomy         # Reformatted taxonomy file using code/python/tax_family.py
    |  +- CDIclear.final.shared             # Final shared file from mothur
    |
    |- code/           # Data analysis scripts
    |  |- R/           # R scipts for figures and analysis
    |      |- day0_rf.R      # Code for random forest analysis for the day of infection
    |      |- day1_rf.R      # Code for random forest analysis for Day 1 of infection
    |      |- preabx_rf.R      # Code for random forest analysis for days prior to antibiotic treatment
    |      +- subsample_shared.R    # Code to identify optimal subsample level for shared file
    |  +- python/       # Python analysis scripts
    |      +- tax_family.py         # Reformats mothur taxonomy files to be more human readable
    |
    |- figures         # All output from workflows and analyses
    |  |- day0_RF_plot.pdf         # Plot for random forest analysis for the day of infection
    |  |- day1_RF_plot.pdf      # Plot for random forest analysis for Day 1 of infection
    |  |- preabx_RF_plot.pdf     # Plotfor random forest analysis for days prior to antibiotic treatment
    |  +- subsample_plot.pdf  # Plot to identify optimal subsample level for shared file


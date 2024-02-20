package ru.biosoft.physicell.sample_projects.virus_macrophage;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;

class Epithelial extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // bookkeeping
        Microenvironment microenvironment = pCell.getMicroenvironment();
        int nVirus = microenvironment.findDensityIndex( "virus" );
        int nInterferon = microenvironment.findDensityIndex( "interferon" );
        int apoptosis_model_index = pCell.phenotype.death.findDeathModelIndex( "Apoptosis" );

        // compare against viral load. Should I commit apoptosis? 
        //        if( pCell.internalizedSubstrates[0] > 0 || pCell.internalizedSubstrates[1] > 0 )
        //        {
        //            double b = 5;
        //        }

        double virus = phenotype.molecular.internalized_total_substrates[nVirus];
        if( virus >= pCell.custom_data.get( "burst_virion_count" ) )
        {
            System.out.println( "\t\tburst!" );
            pCell.lyseCell(); // start_death( apoptosis_model_index );
            pCell.functions.updatePhenotype = null;
            return;
        }

        // replicate virus particles inside me 

        if( virus >= pCell.custom_data.get( "min_virion_count" ) )
        {
            double new_virus = pCell.custom_data.get( "viral_replication_rate" );
            new_virus *= dt;
            phenotype.molecular.internalized_total_substrates[nVirus] += new_virus;
        }

        if( virus >= pCell.custom_data.get( "virion_threshold_for_interferon" ) )
        {
            phenotype.secretion.secretionRates[nInterferon] = pCell.custom_data.get( "max_interferon_secretion_rate" );
        }

        //	 double implicit_Euler_constant = 
        //		(1.0 + dt * pCell.custom_data.get"virus_digestion_rate"] );
        //	phenotype.molecular.internalized_total_substrates[nVirus] /= implicit_Euler_constant; 


        // if I have too many 
        // if I have too many 

    }
}
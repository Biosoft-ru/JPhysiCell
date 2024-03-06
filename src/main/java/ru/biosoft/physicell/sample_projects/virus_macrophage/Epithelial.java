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
        double virus = phenotype.molecular.internSubstrates[nVirus];
        if( virus >= pCell.customData.get( "burst_virion_count" ) )
        {
            pCell.lyseCell(); // start_death( apoptosis_model_index );
            pCell.functions.updatePhenotype = null;
            return;
        }

        // replicate virus particles inside me 
        if( virus >= pCell.customData.get( "min_virion_count" ) )
        {
            double new_virus = pCell.customData.get( "viral_replication_rate" ) * dt;
            phenotype.molecular.internSubstrates[nVirus] += new_virus;
        }

        if( virus >= pCell.customData.get( "virion_threshold_for_interferon" ) )
        {
            phenotype.secretion.secretionRates[nInterferon] = pCell.customData.get( "max_interferon_secretion_rate" );
        }

        //	 double implicit_Euler_constant = 
        //		(1.0 + dt * pCell.custom_data.get"virus_digestion_rate"] );
        //	phenotype.molecular.internalized_total_substrates[nVirus] /= implicit_Euler_constant; 
    }
}
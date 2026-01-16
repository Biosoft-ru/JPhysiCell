package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.micromets.MicrometsModel.LymphNode;

public class dc_history
{

    MicrometsModel model;
    LymphNode lympNode;
    History history;
    
    public dc_history(MicrometsModel model)
    {
        this.model = model;
        this.lympNode = model.lympNode;
        this.history = model.history; 
    }

    void DC_history_model(Cell pCell, Phenotype phenotype, double dt)
    {
        int DC_type = model.getCellDefinition( "DC" ).type;
        // bookkeeping -- find microenvironment variables we need

        // bookkeeping -- find custom data we need
        double DCprob = model.getParameterDouble( "DC_leave_rate" ) * dt;

        // do nothing if dead
        if( phenotype.death.dead == true )
        {
            return;
        }

        // if not DC, do nothing
        if( pCell.type != DC_type )
        {
            return;
        }

        // (Adrianne) if DC is already activated, then check whether it leaves the tissue
        if( pCell.customData.get( "activated_immune_cell" ) > 0.5 && model.getRNG().UniformRandom() < DCprob )
        {
            // (Adrianne) DC leaves the tissue and so we lyse that DC
            System.out.println( "DC leaves tissue" );
            pCell.lyseCell();
            //            #pragma omp critical
            {
                model.lympNode.DCAMOUNT++;
            } // add one
            return;

        }

        return;
    }

    void DC_history_main_model(double dt)
    {
        lympNode.DCAMOUNT = 0;

        //        #pragma omp parallel for
        for( Cell pC : model.getMicroenvironment().getAgents( Cell.class ) )
        {
            //            Cell pC = (*all_cells)[n];
            if( pC.phenotype.death.dead == false )
            {
                DC_history_model( pC, pC.phenotype, dt );
            }
        }
        history.addFront( lympNode.DCAMOUNT );
        //        std::rotate(history.rbegin(),history.rbegin()+1,history.rend());
        //        history.front() = lympNode.DCAMOUNT;

        return;
    }
}

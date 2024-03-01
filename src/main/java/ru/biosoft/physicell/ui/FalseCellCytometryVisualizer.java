package ru.biosoft.physicell.ui;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.PhysiCellConstants;

public class FalseCellCytometryVisualizer extends AgentVisualizer
{
    public Color findColor(Cell pCell)
    {
        //        static std::vector< std::string > output( 4 , "rgb(0,0,0)" );

        // First, check for death. Use standard dead colors and exit


        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.apoptotic ) // Apoptotic - Red
        {
            return new Color( 125, 0, 0 );
            //            output[0] = "rgb(255,0,0)";
            //            output[2] = "rgb(125,0,0)";
            //            return output; 
        }

        // Necrotic - Brown
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
        {
            return new Color( 139, 69, 19 );
            //            output[0] = "rgb(250,138,38)";
            //            output[2] = "rgb(139,69,19)";
            //            return output; 
        }


        // Check if this coloring function even makes sense, and if so,

        if( pCell.phenotype.cycle.code != PhysiCellConstants.flow_cytometry_separated_cycle_model
                && pCell.phenotype.cycle.code != PhysiCellConstants.flow_cytometry_cycle_model )
        {
            return Color.white;
        }

        // G0/G1 and G1 are blue 
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.G0G1_phase
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.G1_phase )
        {
            return new Color( 0, 40, 255 );
            //            output[0] = "rgb(0,80,255)"; 
            //            output[2] = "rgb(0,40,255)"; 
            //            return output; 
        }

        // G0 is pale blue 
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.G0_phase )
        {
            return new Color( 20, 100, 255 );
            //            output[0] = "rgb(40,200,255)";
            //            output[2] = "rgb(20,100,255)";
            //            return output;
        }

        // S is magenta  
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.S_phase )
        {
            return new Color( 190, 0, 190 );
            //            output[0] = "rgb(255, 0, 255)";
            //            output[2] = "rgb(190,0,190)";
            //            return output;
        }

        // G2 is yellow
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.G2_phase )
        {
            return new Color( 190, 190, 0 );
//            output[0] = "rgb(255, 255, 0)";
//            output[2] = "rgb(190, 190, 0)";
//            return output;
        }

        // G2/M and M are green 
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.G2M_phase
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.M_phase )
        {
            return new Color( 0, 190, 0 );
//            output[0] = "rgb(0,255,0)";
//            output[2] = "rgb(0,190,0)";
//
//            return output;
        }

        return Color.white;
    }
}

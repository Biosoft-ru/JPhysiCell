package ru.biosoft.physicell.sample_projects.celltypes3;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.ui.AgentVisualizer;
import ru.biosoft.physicell.xml.ModelReader;

public class RegularAgentVisualizer extends AgentVisualizer
{
    Model model;
    Color aColor;
    Color bColor;
    Color cColor;

    public RegularAgentVisualizer(Model model)
    {
        this.model = model;
        aColor = ModelReader.readColor( model.getParameter( "A_color" ) );
        bColor = ModelReader.readColor( model.getParameter( "B_color" ) );
        cColor = ModelReader.readColor( model.getParameter( "C_color" ) );
    }
    @Override
    public Color findColor(Cell pCell)
    {
        CellDefinition aCD = CellDefinition.getCellDefinition( "A" );
        CellDefinition bCD = CellDefinition.getCellDefinition( "B" );
        CellDefinition cCD = CellDefinition.getCellDefinition( "C" );

        int A_type = aCD.type;
        int B_type = bCD.type;
        int C_type = cCD.type;

        // start with flow cytometry coloring 
        Color output = Color.black;
        //            String[] output = {"black", "black", "black", "black"};

        // color live C 
        if( pCell.type == A_type )
        {
            output = aColor;//ModelReader.readColor( model.getParameter( "A_color" ) );
            //                output[0] = parameters.strings( "A_color" );
            //                output[2] = parameters.strings( "A_color" );
        }

        // color live B
        if( pCell.type == B_type )
        {
            output = bColor;//ModelReader.readColor( model.getParameter( "B_color" ) );
            //                output[0] = parameters.strings( "B_color" );
            //                output[2] = parameters.strings( "B_color" );
        }

        // color live C
        if( pCell.type == C_type )
        {
            output = cColor;//ModelReader.readColor( model.getParameter( "C_color" ) );
            //                output[0] = parameters.strings( "C_color" );
            //                output[2] = parameters.strings( "C_color" );
        }

        if( pCell.phenotype.death.dead == true )
        {
            // Necrotic - Brown
            if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                    || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                    || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
            {
                output = new Color( 123, 63, 0 );//  "chocolate";
            }
            else
            {
                output = Color.black;//"black";
            }
        }
        return output;
    }
}
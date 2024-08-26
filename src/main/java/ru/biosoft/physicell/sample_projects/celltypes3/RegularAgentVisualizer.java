package ru.biosoft.physicell.sample_projects.celltypes3;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
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
        aColor = ModelReader.readColor( model.getParameterString( "A_color" ) );
        bColor = ModelReader.readColor( model.getParameterString( "B_color" ) );
        cColor = ModelReader.readColor( model.getParameterString( "C_color" ) );
    }

    @Override
    public Color[] findColors(Cell cell)
    {
        Color[] output = new Color[] {Color.black, Color.black, Color.black, Color.black};

        if( cell.type == cell.getModel().getCellDefinition( "A" ).type )
        {
            output[0] = aColor;
            output[2] = aColor;
        }
        else if( cell.type == cell.getModel().getCellDefinition( "B" ).type )
        {
            output[0] = bColor;
            output[2] = bColor;
        }
        else if( cell.type == cell.getModel().getCellDefinition( "C" ).type )
        {
            output[0] = cColor;
            output[2] = cColor;
        }

        if( cell.phenotype.death.dead )
        {
            // Necrotic - Brown
            if( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                    || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                    || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
                output[2] = new Color( 123, 63, 0 );//  "chocolate";
            else
                output[2] = Color.black;
        }
        return output;
    }
}
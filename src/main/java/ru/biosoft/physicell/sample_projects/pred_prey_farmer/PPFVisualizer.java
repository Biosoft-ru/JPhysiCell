package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.ui.AgentVisualizer2;

public class PPFVisualizer extends AgentVisualizer2
{
    @Override
    public Color findColor(Cell pCell)
    {
        CellDefinition pFarmerDef = CellDefinition.getCellDefinition( "farmer" );
        CellDefinition pPreyDef = CellDefinition.getCellDefinition( "prey" );
        CellDefinition pPredDef = CellDefinition.getCellDefinition( "predator" );

        if( pCell.type == pFarmerDef.type )
        {
            return Color.gray;
            //{ "grey", "black", "grey", "grey" }; }
        }

        if( pCell.type == pPreyDef.type )
        {
            return Color.blue;
        }//{ "blue", "black", "blue", "blue" }; }

        if( pCell.type == pPredDef.type )
        {
            return Color.orange;
        }
        //            { "orange", "black", "orange", "orange" }; };};

        return super.findBorderColor( pCell );//paint_by_number_cell_coloring( pCell );
    }
}
package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.PhysiCellConstants;

public class AgentVisualizer2 extends AgentVisualizer
{
    List<Color> colors = new ArrayList<>();
    boolean isInit = false;
    public void init()
    {
        colors.add( Color.gray );//"grey" ); // default color will be grey 

        colors.add( Color.red );//"red" );
        colors.add( Color.yellow );///"yellow" );
        colors.add( Color.green );//"green" );
        colors.add( Color.blue );//"blue" );

        colors.add( Color.magenta );//"magenta" );
        colors.add( Color.orange );//"orange" );
        colors.add( new Color( 50, 205, 50 ) );//"lime" );
        colors.add( Color.cyan );//"cyan" );

        colors.add( new Color( 255, 105, 180 ) );//"hotpink" );
        colors.add( new Color( 255, 218, 185 ) );//"peachpuff" );
        colors.add( new Color( 143, 188, 143 ) );//"darkseagreen" );
        colors.add( new Color( 135, 206, 250 ) );//"lightskyblue" );
        isInit = true;
    }

    @Override
    public Color[] findColors(Cell cell)
    {
        if( !isInit )
            init();

        // start all black      
        Color output = Color.black;
        //        std::vector<std::string> output = { "black", "black", "black", "black" }; 

        // paint by number -- by cell type 
        //        std::string interior_color = "white"; 
        Color interiorColor = Color.white;
        if( cell.type < 13 )
        {
            interiorColor = colors.get( cell.type );
        }

        output = interiorColor; // set cytoplasm color     

        // necrotic cells are brown 
        if( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
        {
            interiorColor = new Color( 139, 69, 19 );//"saddlebrown";
        }
        // apoptotic cells are white 
        if( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.apoptotic )
        {
            interiorColor = Color.black;//"black";
        }

        output = interiorColor; // set cytoplasm color 
        //        output[2] = interior_color; // set cytoplasm color 
        //        output[3] = interior_color; // set cytoplasm color 

        //        output[1] = "black";

        return new Color[] {output};
    }
}

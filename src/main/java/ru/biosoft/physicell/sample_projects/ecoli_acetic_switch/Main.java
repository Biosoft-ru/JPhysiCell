package ru.biosoft.physicell.sample_projects.ecoli_acetic_switch;

import java.util.logging.Logger;

import biouml.plugins.fbc.ApacheModel;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.xml.ModelReader;

public class Main
{
    private static String settingsPath = "config/PhysiCell_settings.xml";
    private static String resultPath = "C:/Users/Damag/BIOFVM/projects/ecoli_acetic_switch/result2";
    protected static Logger log = Logger.getLogger( ApacheModel.class.getName() );

    public static void main(String ... strings) throws Exception
    {
        if( strings != null && strings.length > 0 )
            resultPath = strings[0];

        Model model = new ModelReader().read( Main.class.getResourceAsStream( settingsPath ), EcoliAceticSwitch.class );
        double mechanics_voxel_size = 30;
        model.createContainer( mechanics_voxel_size );
        model.setResultFolder( resultPath );
        model.setWriteDensity( true );
        model.setSaveFull( true );
        model.addVisualizer( 0, "oxygen" ).setStubstrateIndex( 0 ).setMaxDensity( 38 );
        model.addVisualizer( 0, "glucose" ).setStubstrateIndex( 1 ).setMaxDensity( 50 );
        model.addVisualizer( 0, "acetate" ).setStubstrateIndex( 2 ).setMaxDensity( 50 );

        model.init();
        System.out.println( model.display() );
        model.simulate();
    }
}
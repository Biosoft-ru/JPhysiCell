package ru.biosoft.biofvm.cell;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;

import ru.biosoft.biofvm.BasicAgent;

public class Output
{
    public static void writePov(Set<Cell> agents, double timepoint, double scale, String fileName) throws IOException
    {
        final int TUMOR_TYPE = 0;
        final int VESSEL_TYPE = 1;
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( fileName ) ) ))
        {
            //        filename = filename.padEnd( 1024, ' ' );
            //            filename = String.format( "%s/cells_%d.pov", PhysiCellСonstants.folder, Math.round( timepoint ) );

            //            FileWriter povFile = new FileWriter( filename );
            bw.write( "#include \"colors.inc\" \n" );
            bw.write( "#include \"header.inc\" \n" );

            for( BasicAgent agent : agents )
            {
                Cell cell = (Cell)agent;
                String _nameCore;
                if( cell.phenotype.cycle != null )
                {
                    int code = cell.phenotype.cycle.currentPhase().code;
                    if( code == PhysiCellConstants.Ki67_positive_premitotic || code == PhysiCellConstants.Ki67_positive_postmitotic
                            || code == PhysiCellConstants.Ki67_positive || code == PhysiCellConstants.Ki67_negative
                            || code == PhysiCellConstants.live )
                        _nameCore = "LIVE";
                    else if( code == PhysiCellConstants.apoptotic )
                        _nameCore = "APOP";
                    else if( code == PhysiCellConstants.necrotic_swelling || code == PhysiCellConstants.necrotic_lysed
                            || code == PhysiCellConstants.necrotic )
                        _nameCore = "NEC";
                    else if( code == PhysiCellConstants.debris )
                        _nameCore = "DEBR";
                    else
                        _nameCore = "MISC";
                }
                else if( cell.type == TUMOR_TYPE )
                    _nameCore = "LIVE";
                else if( cell.type == VESSEL_TYPE )
                    _nameCore = "ENDO";
                else
                    _nameCore = "MISC";
                String center = "<" + Double.toString( cell.position[0] / scale ) + "," + Double.toString( cell.position[1] / scale ) + ","
                        + Double.toString( cell.position[2] / scale ) + ">";
                String core = "sphere {\n\t" + center + "\n\t " + Double.toString( cell.phenotype.geometry.radius / scale )
                        + "\n\t FinishMacro ( " + center + "," + _nameCore + "Finish," + _nameCore + "*1)\n}\n";
                //                povFile.write( core );
                bw.write( core );
            }
            bw.write( "#include \"footer.inc\" \n" );
            //            povFile.write( "#include \"footer.inc\" \n" );
            //            povFile.close();
            //        }
            //        catch( IOException e )
            //        {
            //            e.printStackTrace();
            //        }
            //        return 0;
            //    }
        }
    }

    public static void writeCellReport(Set<Cell> agents, double timepoint, String fileName) throws IOException
    {
        //        std::string filename; 
        //        filename.resize( 1024 ); 
        //  sprintf( (char*) filename.c_str() , "output//cells_%i.txt" , (int)round(timepoint) ); 
        //        sprintf( (char*) filename.c_str() , "%s/cells_%i.txt" , PhysiCell_settings.folder.c_str() , (int)round(timepoint) ); 
        //        std::ofstream povFile (filename.c_str(), std::ofstream::out);
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( fileName ) ) ))
        {
            bw.write(
                    "\tID\tx\ty\tz\tradius\tvolume_total\tvolume_nuclear_fluid\tvolume_nuclear_solid\tvolume_cytoplasmic_fluid\tvolume_cytoplasmic_solid\tvolume_calcified_fraction\tphenotype\telapsed_time\n" );
            int phenotype_code;
            for( BasicAgent agent : agents )
            {
                Cell cell = (Cell)agent;
                phenotype_code = cell.phenotype.cycle.currentPhase().code;
                // phenotype_code = phases.size()>0?cell.phenotype.cycle.phases[cell.phenotype.current_phase_index].code:-1;
                bw.write( cell.ID + "\t" + cell.position[0] + "\t" + cell.position[1] + "\t" + cell.position[2] + "\t" );
                bw.write( +cell.phenotype.geometry.radius + "\t" + cell.phenotype.volume.total + "\t" + cell.phenotype.volume.nuclear_fluid
                        + "\t" + cell.phenotype.volume.nuclear_solid + "\t" + cell.phenotype.volume.cytoplasmic_fluid + "\t"
                        + cell.phenotype.volume.cytoplasmic_solid + "\t" + cell.phenotype.volume.calcified_fraction + "\t" + phenotype_code
                        +
                        // "\t"+ cell.phenotype.cycle.phases[cell.phenotype.current_phase_index].elapsed_time +std::endl;       
                        "\t" + cell.phenotype.cycle.data.elapsedTimePhase );
            }
        }
    }
}

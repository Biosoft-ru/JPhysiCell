package ru.biosoft.biofvm.examples;

import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;

import ru.biosoft.biofvm.BasicAgent;
import ru.biosoft.biofvm.ConstantCoefficientsLOD3D;
import ru.biosoft.biofvm.Microenvironment;
import ru.biosoft.biofvm.SimpleVisualizer;
import ru.biosoft.biofvm.VectorUtil;

/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.0.4) [1]        #
#                                                                           #
# [1] A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient  #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics, 2015. DOI: 10.1093/bioinformatics/btv730               #
#############################################################################
#                                                                           #
# Copyright 2015 Paul Macklin and the BioFVM Project                        #
#                                                                           #
# Licensed under the Apache License, Version 2.0 (the "License");           #
# you may not use this file except in compliance with the License.          #
# You may obtain a copy of the License at                                   #
#                                                                           #
#    http://www.apache.org/licenses/LICENSE-2.0                             #
#                                                                           #
# Unless required by applicable law or agreed to in writing, software       #
# distributed under the License is distributed on an "AS IS" BASIS,         #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #
# See the License for the specific language governing permissions and       #
# limitations under the License.                                            #
#############################################################################
*/
public class Tutorial1
{
    //Grid options
    private static int size = 1000; // x size = y size = z size
    private static int cellSize = 10; // x size = y size = z size

    //Substrate options
    private static int diffusionCoefficients = 1000;// microns^2 / min 
    private static double decayRate = 0.01;
    private static double saturationDensity = 1;

    //Agents options
    private static int sourceAgentsNumber = 500;
    private static int sinkAgentsNumber = 500;
    private static double agentRadius = 5;
    private static boolean generateAtSlice = true; //create agents only at slice we watch
	private static int agentSecretionRate = 10;
	private static double agentUptakeRate = 0.8;

    //Simulation options
    private static double tMax = 10;
    private static double dt = 0.01;
    private static String resultPath = "C:/Users/Damag/BIOFVM/";
    private static int zSlice = 500; //watch slice in the middle (for images)
    private static boolean outputImages = true;
    private static boolean outputGIF = true;
    private static int imgTimeDelta = 10; // time delta between images
    private static boolean outputTables = false;

    private static ImageWriter writer;
    private static ImageOutputStream ios;

    public static void main(String ... args) throws Exception
    {
        double timeStart = System.currentTimeMillis();

        Microenvironment m = new Microenvironment();
        m.name = "substrate scale";
        m.set_density( 0, "substrate1", "dimensionless" );
        m.spatial_units = "microns";
        m.mesh.units = "microns";
        m.time_units = "minutes";

        m.resize_space_uniform( 0, size, 0, size, 0, size, cellSize );
        m.displayInformation();

        //set initial density - more in the middle, gradually decrease
        double[] center = new double[3];
        center[0] = ( m.mesh.bounding_box[0] + m.mesh.bounding_box[3] ) / 2;
        center[1] = ( m.mesh.bounding_box[1] + m.mesh.bounding_box[4] ) / 2;
        center[2] = ( m.mesh.bounding_box[2] + m.mesh.bounding_box[5] ) / 2;
        double stddevSquared = size / cellSize;
        stddevSquared *= -stddevSquared;
        for( int i = 0; i < m.number_of_voxels(); i++ )
        {
            double[] displacement = VectorUtil.newDiff( m.voxels( i ).center, center );
            double coeff = VectorUtil.norm_squared( displacement ) / stddevSquared;
            m.density_vector( i )[0] = Math.exp( coeff );
        }

        // register the diffusion solver    
        m.setSolver( new ConstantCoefficientsLOD3D() );

        // register substrates properties 
        m.diffusion_coefficients[0] = diffusionCoefficients; // microns^2 / min 
        m.decay_rates[0] = decayRate;
        m.options.track_internalized_substrates_in_each_agent = false;

        //populate model
        for( int i = 0; i < sourceAgentsNumber; i++ )
        {
            //substrate source agents
            double[] pos = VectorUtil.random( 3, 0, size );
            if( generateAtSlice )
                pos[2] = zSlice; //keep agents at z slice

            BasicAgent agentSource = BasicAgent.createBasicAgent();
            agentSource.registerMicroenvironment( m );
            agentSource.assignPosition( pos );
            agentSource.setRadius( agentRadius );
            agentSource.secretionRates[0] = agentSecretionRate;
            agentSource.saturationDensities[0] = saturationDensity;
            agentSource.setUptakeConstants( dt );
            agentSource.setTag( "Source" );
        }

        for( int i = 0; i < sinkAgentsNumber; i++ )
        {
            double[] pos = VectorUtil.random( 3, 0, size );
            if( generateAtSlice )
                pos[2] = zSlice; //keep agents at z slice 

            BasicAgent agentSink = BasicAgent.createBasicAgent();
            agentSink.registerMicroenvironment( m );
            agentSink.assignPosition( pos );
            agentSink.setRadius( agentRadius );
            agentSink.uptakeRates[0] = agentUptakeRate;
            agentSink.setUptakeConstants( dt );
            agentSink.setTag( "Sink" );
        }
        //visualization of all layers
        //        for( int z = 0; z < 100; z++ )
        //            SimpleVisualizer.draw( m, z, 1000, 1000, resultPath + "/init_slice" + z + ".png" );

        if( outputGIF )
        {
            writer = ImageIO.getImageWritersByFormatName( "GIF" ).next();
            ios = ImageIO.createImageOutputStream( new File( resultPath + "/result.gif" ) );
            writer.setOutput( ios );
            writer.prepareWriteSequence( null );
        }

        if( outputTables )
            m.write_to_matlab( resultPath + "/initial.txt" );
        if( outputImages )
        {
            BufferedImage nextImg = SimpleVisualizer.draw( m, zSlice, 0, resultPath + "/initial.png" );
            if( outputGIF )
                writer.writeToSequence( new IIOImage( nextImg, null, null ), writer.getDefaultWriteParam() );
        }
        double t = 0.0;
        int counter = 0;
        //Simulation starts
        double timeSimulationStart = System.currentTimeMillis();
        while( t < tMax )
        {
            m.simulate_cell_sources_and_sinks( dt );
            m.simulate_diffusion_decay( dt );
            t += dt;
            counter++;

            if( outputTables )
                m.write_to_matlab( resultPath + "/step_" + counter + ".txt" );

            if( outputImages && counter % imgTimeDelta == 0 )
            {
                BufferedImage nextImg = SimpleVisualizer.draw( m, zSlice, Math.round( t * 100 ) / 100.0,
                        resultPath + "/step_" + counter + ".png" );
                if( outputGIF )
                    writer.writeToSequence( new IIOImage( nextImg, null, null ), writer.getDefaultWriteParam() );
            }
        }

        if( outputGIF )
        {
            writer.endWriteSequence();
            writer.reset();
            writer.dispose();
            ios.flush();
            ios.close();
        }
        double totalTime = ( System.currentTimeMillis() - timeStart ) / 1000;
        double simulationTime = ( System.currentTimeMillis() - timeSimulationStart ) / 1000;
        System.out.println( "Done! Elapsed time: " + totalTime + " s. Simulation time: " + simulationTime + " s." );
    }
}

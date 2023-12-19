package ru.biosoft.biofvm.cell;

import java.util.ArrayList;
import java.util.List;

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
public class Death implements Cloneable
{
    public List<Double> rates;
    List<CycleModel> models;
    List<DeathParameters> parameters;

    boolean dead;
    int current_death_model_index;

    public Death()
    {
        rates = new ArrayList<>();// double[0];//.resize( 0 ); 
        models = new ArrayList<>();//= new CycleModel[0];//.resize( 0 ); 
        parameters = new ArrayList<>();//= new DeathParameters[0];//.resize( 0 ); 

        dead = false;
        current_death_model_index = 0;
    }

    int add_death_model(double rate, CycleModel pModel)
    {
        rates.add( rate );
        models.add( pModel );
        parameters.add( new DeathParameters() );
        //    parameters.resize( rates.size() ); //TODO:
        return rates.size() - 1;
    }

    int add_death_model(double rate, CycleModel pModel, DeathParameters deathParameters)
    {
        //    rates.push_back( rate );
        //    models.push_back( pModel ); 
        //    parameters.push_back( death_parameters ); 
        rates.add( rate );
        models.add( pModel );
        parameters.add( deathParameters );
        return rates.size() - 1;
    }

    int find_death_model_index(int code)
    {
        for( int i = 0; i < models.size(); i++ )
        {
            if( models.get( i ).code == code )
            {
                return i;
            }
        }
        return 0; //TODO??
    }

    public int find_death_model_index(String name)
    {
        for( int i = 0; i < models.size(); i++ )
        {
            if( models.get( i ).name == name )
            {
                return i;
            }
        }
        return 0;
    }

    boolean check_for_death(double dt)
    {
        // If the cell is already dead, exit. 
        if( dead == true )
        {
            return false;
        }

        // If the cell is alive, evaluate all the 
        // death rates for each registered death type. 
        int i = 0;
        while( !dead && i < rates.size() )
        {
            double rate = rates.get( i );
            if( rate != 0 && Math.random() < rate * dt )
            {
                // update the Death data structure 
                dead = true;
                current_death_model_index = i;

                // and set the cycle model to this death model 

                return dead;
            }
            i++;
        }

        return dead;
    }

    public void trigger_death(int death_model_index)
    {
        dead = true;
        current_death_model_index = death_model_index;

        /*  
            // if so, change the cycle model to the current death model 
            phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() ); 
            
            // also, turn off motility.
            
            phenotype.motility.is_motile = false; 
            phenotype.motility.motility_vector.assign( 3, 0.0 ); 
            functions.update_migration_bias = NULL;
            
            // turn off secretion, and reduce uptake by a factor of 10 
            phenotype.secretion.set_all_secretion_to_zero();
            phenotype.secretion.scale_all_uptake_by_factor( 0.10 );
            
            // make sure to run the death entry function 
            if( phenotype.cycle.current_phase().entry_function )
            {
        phenotype.cycle.current_phase().entry_function( this, phenotype, dt_ ); 
            }
        */
    }

    public DeathParameters current_parameters()
    {
        return parameters.get( current_death_model_index );
    }

    public CycleModel current_model()
    {
        return models.get( current_death_model_index );
    }

    //double ref
    //    double apoptosis_rate()
    //    {
    //        int nApoptosis = find_death_model_index( PhysiCellConstants.apoptosis_death_model );
    //        return rates[nApoptosis];
    //    }
    //
    //    //double ref
    //    double necrosis_rate()
    //    {
    //        int nNecrosis = find_death_model_index( PhysiCellConstants.necrosis_death_model );
    //        return rates[nNecrosis];
    //    }

    @Override
    public Death clone()
    {
        try
        {
            Death result = (Death)super.clone();
            result.rates = new ArrayList<>( rates );//this.rates  public List<Double> rates; 
            result.parameters = new ArrayList<>();
            for( int i = 0; i < parameters.size(); i++ )
                        result.parameters.add( parameters.get( i ).clone() );
                    result.models = new ArrayList<>();
            for( int i = 0; i < models.size(); i++ )
                result.models.add( models.get( i ).clone() );
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}

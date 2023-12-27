package ru.biosoft.biofvm.cell;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
public class CycleModel implements Cloneable
{

    /* this maps the end_phase_index to the link index in each phase_links[i]
    * So, index_inverse_map[i][j] = k, corresponds to phases[i], phase_links[i][k] (which links from phase i to phase j)
    * transition_rates[i][k] (the transition rate from phase i to phase j
    */
    List<Map<Integer, Integer>> inverse_index_maps;
    String name;
    int code;
    List<Phase> phases;
    List<List<PhaseLink>> phase_links;
    int default_phase_index;
    public CycleData data; // this will be copied to individual cell agents 

    public CycleModel()
    {
        inverse_index_maps = new ArrayList<>();//HashMap<Integer, Integer>();//.resize( 0 );
        name = "unnamed";
        phases = new ArrayList<>();
        phase_links = new ArrayList<>();
        data = new CycleData( this );// = this;
        code = PhysiCellConstants.custom_cycle_model;
        default_phase_index = 0;
    }

    public Phase currentPhase()
    {
        return phases.get( data.currentPhaseIndex );
    }

    public int add_phase(int code, String name)
    {
        int n = phases.size();
        phases.add( new Phase( name, n, code ) );
        phase_links.add( new ArrayList<PhaseLink>() );
        inverse_index_maps.add( new HashMap<Integer, Integer>() );
        data.sync_to_cycle_model();
        return n;
    }

    public int add_phase_link(int start_index, int end_index, PhaseArrest arrestFunction) throws Exception
    {
        List<PhaseLink> links = phase_links.get( start_index );
        int n = links.size();
        PhaseArrest arrest = arrestFunction == null ? null : arrestFunction.getClass().newInstance();
        Phase start = phases.get( start_index );
        Phase end = phases.get( end_index );
        links.add( new PhaseLink( start, start_index, end, end_index, arrest ) );
        inverse_index_maps.get( start_index ).put( end_index, n );//link from start_index to end_index have index n        
        data.sync_to_cycle_model(); // lastly, make sure the transition rates are the right size;
        return n;
    }

    public int add_phase_link(int start_index, int end_index, double rate, PhaseArrest arrestFunction) throws Exception
    {
        int n = add_phase_link( start_index, end_index, arrestFunction );
        data.setTransitionRate( start_index, end_index, rate );
        return n;
    }

    public int findPhaseIndex(int code)
    {
        for( int i = 0; i < phases.size(); i++ )
        {
            if( phases.get( i ).code == code )
            {
                return i;
            }
        }
        return -1;//?
    }

    int find_phase_index(String name)
    {
        for( int i = 0; i < phases.size(); i++ )
        {
            if( phases.get( i ).name == name )
            {
                return i;
            }
        }
        return -1;
    }

    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Cycle Model: " + name + ". PhyCell code: " + code + "\n" );
        sb.append( "Phases and links: (* denotes phase with cell division)\n" );
        for (int i=0; i<phases.size(); i++)
        {
            sb.append( "Phase " + i + " (" + phases.get( i ).name + ")" );
            if( phases.get( i ).divisionAtExit )
                sb.append( "*" );
            sb.append( "\n" );
            sb.append( " links to: \n" );

            for( int k = 0; k < phase_links.get( i ).size(); k++ )
            {
                int j = phase_links.get( i ).get( k ).endPhaseIndex;
                sb.append( "\tPhase " + " (" + phases.get( j ).name + "( with rate " + data.getTransitionRate( i, j ) + " "
                        + data.timeUnits + "^-1; \n" );
            }
            sb.append( "\n" );
        }
        sb.append( "\n" );
        return sb.toString();
    }


    //    double& ::transition_rate( int start_index , int end_index )
    public double transition_rate(int start_index, int end_index)
    {
        return data.getTransitionRate( start_index, end_index );
    }

    public void setTransitionRate(int start_index, int end_index, double rate)
    {
        data.setTransitionRate( start_index, end_index, rate );
    }

    public PhaseLink phase_link(int start_index, int end_index)
    {
        return phase_links.get( start_index ).get( inverse_index_maps.get( start_index ).get( end_index ) );
    }

    public void advance(Cell pCell, Phenotype phenotype, double dt)
    {
        int i = phenotype.cycle.data.currentPhaseIndex;

        phenotype.cycle.data.elapsedTimePhase += dt;

        // Evaluate each linked phase: advance to that phase IF probabiltiy is in the range, and if the arrest function (if any) is false 
        List<PhaseLink> links = phase_links.get( i );
        int j;
        for( int k = 0; k < links.size(); k++ )
        {
            PhaseLink link = links.get( k );
            j = link.endPhaseIndex;
            double transition_rate = phenotype.cycle.data.transitionRates.get( i ).get( k );

            // check for arrest. If arrested, skip to the next transition
            boolean transition_arrested = false;
            if( link.arrestFunction != null )
            {
                transition_arrested = link.arrestFunction.isArrested( pCell, phenotype, dt );
            }
            if( !transition_arrested )
            {
                // check to see if we should transition 
                boolean continue_transition = false;
                if( link.fixedDuration )
                {
                    if( phenotype.cycle.data.elapsedTimePhase > 1.0 / transition_rate )
                    {
                        continue_transition = true;
                    }
                }
                else
                {
                    double prob = transition_rate * dt;
                    if( Math.random() < prob )
                    {
                        continue_transition = true;
                    }
                }

                // if we should transition, check if we're not supposed to divide or die 
                if( continue_transition )
                {
                    // if the phase transition has an exit function, execute it
                    if( link.exitFunction != null )
                    {
                        link.exitFunction.execute( pCell, phenotype, dt );
                    }

                    // check if division or removal are required 
                    if( phases.get( i ).divisionAtExit )
                    {
                        phenotype.flagged_for_division = true;
                    }
                    if( phases.get( i ).removalAtExit )
                    {
                        phenotype.flagged_for_removal = true;
                        return;
                    }
                    // move to the next phase, and reset the elapsed time 
                    phenotype.cycle.data.currentPhaseIndex = j;
                    phenotype.cycle.data.elapsedTimePhase = 0.0;

                    // if the new phase has an entry function, execute it 
                    if( phases.get( j ).entryFunction != null )
                    {
                        phases.get( j ).entryFunction.execute( pCell, phenotype, dt );
                    }
                    return;
                }
            }
        }
    }

    @Override
    public CycleModel clone()
    {
        CycleModel result = new CycleModel();
        result.name = name;
        result.code = code;
        result.phases = new ArrayList<Phase>( phases );
        result.phase_links = this.phase_links; //TODO: check
        result.default_phase_index = this.default_phase_index;
        result.inverse_index_maps = this.inverse_index_maps;//TODO: check
        result.data = this.data.clone( result );
        return result;
    }
}
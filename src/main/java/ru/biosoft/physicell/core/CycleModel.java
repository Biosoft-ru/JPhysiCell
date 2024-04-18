package ru.biosoft.physicell.core;

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
    List<Map<Integer, Integer>> startToEndToLink;
    public String name;
    public int code;
    public List<Phase> phases;
    public List<List<PhaseLink>> phaseLinks;
    int default_phase_index;
    public CycleData data; // this will be copied to individual cell agents 

    public CycleModel()
    {
        startToEndToLink = new ArrayList<>();
        name = "unnamed";
        phases = new ArrayList<>();
        phaseLinks = new ArrayList<>();
        data = new CycleData( this );
        code = PhysiCellConstants.custom_cycle_model;
        default_phase_index = 0;
    }

    public Phase currentPhase()
    {
        return phases.get( data.currentPhaseIndex );
    }

    public int addPhase(int code, String name)
    {
        int n = phases.size();
        phases.add( new Phase( name, n, code ) );
        phaseLinks.add( new ArrayList<PhaseLink>() );
        startToEndToLink.add( new HashMap<Integer, Integer>() );
        data.syncToCycle();
        return n;
    }

    public int addPhaseLink(int startIndex, int endIndex, PhaseArrest arrestFunction) throws Exception
    {
        List<PhaseLink> links = phaseLinks.get( startIndex );
        int n = links.size();
        PhaseArrest arrest = arrestFunction == null ? null : arrestFunction.getClass().newInstance();
        Phase start = phases.get( startIndex );
        Phase end = phases.get( endIndex );
        links.add( new PhaseLink( start, startIndex, end, endIndex, arrest ) );
        startToEndToLink.get( startIndex ).put( endIndex, n );//link from start_index to end_index have index n        
        data.syncToCycle(); // lastly, make sure the transition rates are the right size;
        return n;
    }

    public int addPhaseLink(int startIndex, int endIndex, double rate, PhaseArrest arrestFunction) throws Exception
    {
        int n = addPhaseLink( startIndex, endIndex, arrestFunction );
        data.setBasicTransitionRate( startIndex, endIndex, rate );
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

    int findPhaseIndex(String name)
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
        return display();
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Cycle Model: " + name + " (" + code + ")\n" );
        sb.append( "--------------------------------" );

        if( phases.size() == 1 )
        {
            String name = phases.get( 0 ).name;
            if( phases.get( 0 ).divisionAtExit )
                name += "*";
            sb.append( "\n\t" + name );
            if( phaseLinks.get( 0 ).size() == 1 )
            {
                double duration = PhysiCellUtilities.print( 1.0 / data.getTransitionRate( 0, 0 ) );
                sb.append( " -> " + name + ", duration: " + duration + " " + data.timeUnits );
            }
            return sb.toString();
        }

        for( int i = 0; i < phases.size(); i++ )
        {
            sb.append( "\n\t" + phases.get( i ).name + " " );
            if( phases.get( i ).divisionAtExit )
                sb.append( "*" );

            for( int k = 0; k < phaseLinks.get( i ).size(); k++ )
            {
                int j = phaseLinks.get( i ).get( k ).endPhaseIndex;
                double duration = PhysiCellUtilities.print( 1.0 / data.getTransitionRate( i, j ) );
                if( phaseLinks.get( i ).size() > 1 )
                    sb.append( "\n" );
                sb.append( "-> " + phases.get( j ).name );
                if( phases.get( j ).divisionAtExit )
                    sb.append( "*" );
                sb.append( ", duration: " + duration + " " + data.timeUnits );
            }
        }
        return sb.toString();
    }

    public double transition_rate(int start_index, int end_index)
    {
        return data.getTransitionRate( start_index, end_index );
    }

    /**
     * Sets basic transition rate 
     */
    public void setBasicTransitionRate(int start_index, int end_index, double rate)
    {
        data.setBasicTransitionRate( start_index, end_index, rate );
        data.setTransitionRate( start_index, end_index, rate );
    }

    public void setTransitionRate(int start_index, int end_index, double rate)
    {
        data.setTransitionRate( start_index, end_index, rate );
    }

    public PhaseLink phase_link(int start_index, int end_index)
    {
        return phaseLinks.get( start_index ).get( startToEndToLink.get( start_index ).get( end_index ) );
    }

    public void advance(Cell pCell, Phenotype phenotype, double dt)
    {
        int i = phenotype.cycle.data.currentPhaseIndex;

        phenotype.cycle.data.elapsedTimePhase += dt;

        // Evaluate each linked phase: advance to that phase IF probabiltiy is in the range, and if the arrest function (if any) is false 
        List<PhaseLink> links = phaseLinks.get( i );
        int j;
        for( int k = 0; k < links.size(); k++ )
        {
            PhaseLink link = links.get( k );
            j = link.endPhaseIndex;
            double transitionRate = phenotype.cycle.data.transitionRates.get( i ).get( k );

            // check for arrest. If arrested, skip to the next transition
            boolean transitionArrested = false;
            if( link.arrestFunction != null )
            {
                transitionArrested = link.arrestFunction.isArrested( pCell, phenotype, dt );
            }
            if( !transitionArrested )
            {
                // check to see if we should transition 
                boolean continueTransition = false;
                if( link.fixedDuration )
                {
                    if( phenotype.cycle.data.elapsedTimePhase > 1.0 / transitionRate )
                    {
                        continueTransition = true;
                    }
                }
                else
                {
                    double prob = transitionRate * dt;
//                    double r = PhysiCellUtilities.UniformRandom();
                    //                    if( pCell.type_name.equals( "A" ) )
                    //                        System.out.println( pCell + " " + prob + " " + r );
                    if( PhysiCellUtilities.checkRandom(prob))// r < prob )
                    {
                        //                        System.out.println( "Random " + r + " , probability " + prob );
                        continueTransition = true;
                    }
                }

                // if we should transition, check if we're not supposed to divide or die 
                if( continueTransition )
                {
                    // if the phase transition has an exit function, execute it
                    if( link.exitFunction != null )
                    {
                        link.exitFunction.execute( pCell, phenotype, dt );
                    }

                    // check if division or removal are required 
                    if( phases.get( i ).divisionAtExit )
                    {
                        phenotype.flaggedForDivision = true;
                    }
                    if( phases.get( i ).removalAtExit )
                    {
                        phenotype.flaggedForRemoval = true;
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
        result.phaseLinks = this.phaseLinks; //TODO: check
        result.default_phase_index = this.default_phase_index;
        result.startToEndToLink = this.startToEndToLink;//TODO: check
        result.data = this.data.clone( result );
        return result;
    }
}
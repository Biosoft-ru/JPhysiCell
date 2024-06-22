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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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
public class HypothesisRule
{
    Map<String, Integer> signalsMap = new HashMap<>();

    String cellType;
    CellDefinition cd;

    String behavior;
    double baseValue;
    double maxValue;
    double minValue;

    List<String> signals;
    List<Boolean> responses;
    List<Double> halfMaxes;
    List<Double> hillPowers;
    List<Boolean> appliesToDead;

    List<String> upSignals;
    List<Double> upHalfMaxes;
    List<Double> upHillPowers;
    List<Boolean> upAppliesToDead;

    List<String> downSignals;
    List<Double> downHalfMaxes;
    List<Double> downHillPowers;
    List<Boolean> downAppliesToDead;

    HypothesisRule()
    {
        signalsMap.clear();

        behavior = "none";
        baseValue = 1.0;
        maxValue = 10.0;
        minValue = 0.1;

        signals = new ArrayList<>();
        responses = new ArrayList<>();
        halfMaxes = new ArrayList<>();
        hillPowers = new ArrayList<>();
        appliesToDead = new ArrayList<>();

        upSignals = new ArrayList<>();
        upHalfMaxes = new ArrayList<>();
        upHillPowers = new ArrayList<>();
        upAppliesToDead = new ArrayList<>();

        downSignals = new ArrayList<>();
        downHalfMaxes = new ArrayList<>();
        downHillPowers = new ArrayList<>();
        downAppliesToDead = new ArrayList<>();

        cellType = "none";
        cd = null;
    }

    String convert(boolean input)
    {
        return ( input ) ? "increases" : "decreases";
    }

    String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "For cell type " + cellType + ": \n" );
        for( int j = 0; j < signals.size(); j++ )
            sb.append( "\n" + signals.get( j ) + " " + convert( responses.get( j ) ) + " " + behavior );
        return sb.toString();
    }

    String reducedDisplay()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
            sb.append( "\n" + signals.get( j ) + " " + convert( responses.get( j ) ) + " " + behavior );
        return sb.toString();
    }

    String detailedDisplay()
    {
        StringBuilder sb = new StringBuilder();
        // sb.append( "For cell type " + cell_type + ": " + std::endl; 
        sb.append( behavior + " is modulated from " + minValue + " to " + maxValue + " with a base value of " + baseValue + "\n" );
        sb.append( "--------------------------------------------------------\n" );
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( "\t" + signals.get( j ) + " " + convert( responses.get( j ) ) + " " + behavior + " with half-max "
                    + halfMaxes.get( j ) + " and Hill power " + hillPowers.get( j ) + "." );
            if( appliesToDead.get( j ) )
            {
                sb.append( " Rule applies to dead cells." );
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }

    String English_detailed_display()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( signals.get( j ) + " " + convert(responses.get( j )) );
            sb.append( behavior + " from " + baseValue + " towards " );
            if( responses.get( j ) )
                sb.append( maxValue );
            else
                sb.append( minValue );
            sb.append( " with a Hill response, with half-max " + halfMaxes.get( j ) );
            sb.append( " and Hill power " + hillPowers.get( j ) + "." );
            if( appliesToDead.get( j ) )
                sb.append( " Rule applies to dead cells." );
            sb.append( "\n" );
        }
        return sb.toString();
    }

    String English_detailed_display_HTML()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( "<li>" + signals.get( j ) + " " + convert( responses.get( j ) ) );
            sb.append( behavior + " from " + baseValue + " towards " );
            if( responses.get( j ) )
                sb.append( maxValue );
            else
                sb.append( minValue );
            sb.append( " with a Hill response, with half-max " + halfMaxes.get( j ) );
            sb.append( " and Hill power " + hillPowers.get( j ) + "." );
            if( appliesToDead.get( j ) )
                sb.append( " Rule applies to dead cells." );
            sb.append( "</li>\n" );
        }
        return sb.toString();
    }

    String English_display()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( signals.get( j ) + " " + convert( responses.get( j ) ) );
            sb.append( behavior + "\n" );
        }
        return sb.toString();
    }

    String English_display_HTML()
    {
        StringBuilder sb = new StringBuilder();
        for( int j = 0; j < signals.size(); j++ )
        {
            sb.append( "<li>" + signals.get( j ) + " " + convert( responses.get( j ) ) );
            sb.append( behavior + "</li>\n" );
        }
        return sb.toString();
    }

    void addSignal(String signal, double half_max, double hill_power, String response) throws Exception
    {
        // check: is this a valid signal? (is it in the dictionary?)
        if( SignalBehavior.findSignalIndex( signal ) < 0 )
        {
            throw new Exception( "Warning! Attempted to add signal " + signal + " which is not in the dictionary."
                    + "Either fix your model or add the missing signal to the simulation." );
        }

        // check to see if it's already there 
        int n = findSignal( signal );

        // if so, then just warn and exit.  
        if( n > -1 )
        {
            System.out.println( "Warning! Signal " + signal + " was already part of the rule. Ignoring input." );
            return;
        }

        // add the signal; 
        signalsMap.put( signal, signalsMap.size() );

        boolean bResponse = false; // true if up-regulate, false if down
        if( response.equals( "increase" ) || response.equals( "increases" ) || response.equals( "promotes" ) )
            bResponse = true;

        signals.add( signal );
        halfMaxes.add( half_max );
        hillPowers.add( hill_power );
        responses.add( bResponse );
        appliesToDead.add( false );

        // separate into up and down for our convenience 
        if( bResponse )
        {
            upSignals.add( signal );
            upHalfMaxes.add( half_max );
            upHillPowers.add( hill_power );
            upAppliesToDead.add( false );
        }
        else
        {
            downSignals.add( signal );
            downHalfMaxes.add( half_max );
            downHillPowers.add( hill_power );
            downAppliesToDead.add( false );
        }
    }

    void addSignal(String signal, String response) throws Exception
    {
        addSignal( signal, 0.5, 3.0, response );
    }

    double evaluate(double[] signalValues, boolean dead)
    {
        List<Double> upSignal = new ArrayList<>();
        List<Double> downSignal = new ArrayList<>();
        boolean applyRule = false;

        // need to modify to evaluate if cell is live or if the rule is allowed for dead cells 
        for( int j = 0; j < signalValues.length; j++ )
        {
            if( appliesToDead.get( j ) || !dead )
            {
                if( responses.get( j ) )
                    upSignal.add( signalValues[j] );
                else
                    downSignal.add( signalValues[j] );
                applyRule = true;
            }
            else
            {
                // new oin sep 7 , 2022
                if( responses.get( j ) )
                    upSignal.add( 0.0 );
                else
                    downSignal.add( 0.0 );
            }
        }

        // March 27, 2023
        // if none of the rules apply, the return an absurdly low value 
        // to signal that the parameter value shoudl not be written 
        if( !applyRule )
            return -9e99;

        // up-regulation part 
        double hu = BasicSignaling.multivariate_Hill_response_function( upSignal, upHalfMaxes, upHillPowers );
        double u = baseValue + ( maxValue - baseValue ) * hu;

        // then the down-regulation part 
        double du = BasicSignaling.multivariate_Hill_response_function( downSignal, downHalfMaxes, downHillPowers );
        double output = u + ( minValue - u ) * du;

        return output;
    }

    double evaluate(double[] values)
    {
        return evaluate( values, true );
    }

    double evaluate(Cell cell) throws Exception
    {
        if( !cell.phenotype.cycle.name.equals( "Flow cytometry model (separated)" ) )
        {
            double a = 5;
        }
        double[] signalValues = new double[signals.size()];
        for( int i = 0; i < signals.size(); i++ )
            signalValues[i] = SignalBehavior.getSingleSignal( cell, signals.get( i ) );

        // now, get live/dead value 
        boolean dead = SignalBehavior.getSingleSignal( cell, "dead" ) != 0;
        double out = evaluate( signalValues, dead );

        // new March 27, 2023 
        // if the rule was found to not apply, then just get the prior value 
        if( out < -9e90 )
            out = SignalBehavior.getSinglBehavior( cell, this.behavior );
        return out;
    }

    void apply(Cell cell) throws Exception
    {
        double param = evaluate( cell );
        SignalBehavior.setSingleBehavior( cell, behavior, param );
    }

    void sync(Model model, CellDefinition pCD)
    {
        if( pCD == null )
            return;
        cellType = pCD.name;
        SignalBehavior.getSingleBaseBehavior( model, pCD, behavior );
    }

    void sync(Model model, String cellName)
    {
        sync( model, model.getCellDefinition( cellName ) );
    }

    int findSignal(String name)
    {
        if( signalsMap.containsKey( name ) )
            return signalsMap.get( name );
        return -1;
    }

    void setHalfMax(String name, double hm)
    {
        int n = findSignal( name );
        if( n < 0 )
            return;

        halfMaxes.set( n, hm );

        if( responses.get( n ) )
        {
            for( int m = 0; m < upSignals.size(); m++ )
            {
                if( upSignals.get( m ).equals( name ) )
                {
                    upHalfMaxes.set( m, hm );
                }
            }
        }
        else
        {
            for( int m = 0; m < downSignals.size(); m++ )
            {
                if( downSignals.get( m ).equals( name ) )
                {
                    downHalfMaxes.set( m, hm );
                }
            }
        }
    }

    void setHillPower(String name, double hp)
    {
        int n = findSignal( name );
        if( n < 0 )
            return;

        hillPowers.set( n, hp );
        if( responses.get( n ) )
        {
            for( int m = 0; m < upSignals.size(); m++ )
            {
                if( upSignals.get( m ).equals( name ) )
                {
                    upHillPowers.set( m, hp );
                }
            }
        }
        else
        {
            for( int m = 0; m < downSignals.size(); m++ )
            {
                if( downSignals.get( m ).equals( name ) )
                {
                    downHillPowers.set( m, hp );
                }
            }
        }
    }

    void setResponse(String name, String response)
    {
        int n = findSignal( name );
        if( n < 0 )
            return;

        boolean bResponse = false; // true if up-regulate, false if down
        if( response.equals( "increase" ) || response.equals( "increases" ) || response.equals( "promotes" ) )
            bResponse = true;

        // this is already my response? if so exit 
        if( bResponse == responses.get( n ) )
            return;

        if( responses.get( n ) )//[n] == true )
        {
            // need to switch from up to down 
            // find current index 
            int ci = -1;
            for( int m = 0; m < upSignals.size(); m++ )
            {
                if( upSignals.get( m ).equals( name ) )//[m] == name )
                {
                    ci = m;
                }
            }
            // swap last inot that position 
            upHalfMaxes.set( ci, upHalfMaxes.get( upHalfMaxes.size() - 1 ) );
            upHillPowers.set( ci, upHillPowers.get( upHillPowers.size() - 1 ) );
            upSignals.set( ci, upSignals.get( upSignals.size() - 1 ) );
            upAppliesToDead.set( ci, upAppliesToDead.get( upAppliesToDead.size() - 1 ) );

            // reduce size by one
            upHalfMaxes.remove( upHalfMaxes.size() - 1 );
            upHillPowers.remove( upHillPowers.size() - 1 );
            upSignals.remove( upSignals.size() - 1 );
            upAppliesToDead.remove( upAppliesToDead.size() - 1 );

            // move to the other side 
            downHalfMaxes.add( halfMaxes.get( n ) );
            downHillPowers.add( hillPowers.get( n ) );
            downSignals.add( signals.get( n ) );
            downAppliesToDead.add( appliesToDead.get( n ) );
        }
        else
        {
            // need to switch from down to up find current index 
            int ci = -1;
            for( int m = 0; m < downSignals.size(); m++ )
            {
                if( downSignals.get( m ).equals( name ) )
                {
                    ci = m;
                }
            }
            // swap last inot that position 
            downHalfMaxes.set( ci, downHalfMaxes.get( downHalfMaxes.size() - 1 ) );
            downHillPowers.set( ci, downHillPowers.get( downHillPowers.size() - 1 ) );
            downSignals.set( ci, downSignals.get( downSignals.size() - 1 ) );
            downAppliesToDead.set( ci, downAppliesToDead.get( downAppliesToDead.size() - 1 ) );

            // reduce size by one
            downHalfMaxes.remove( downHalfMaxes.size() - 1 );
            downHillPowers.remove( downHillPowers.size() - 1 );
            downSignals.remove( downSignals.size() - 1 );
            downAppliesToDead.remove( downAppliesToDead.size() - 1 );

            // move to the other side 
            upHalfMaxes.add( halfMaxes.get( n ) );
            upHillPowers.add( hillPowers.get( n ) );
            upSignals.add( signals.get( n ) );
            upAppliesToDead.add( appliesToDead.get( n ) );
        }
        responses.set( n, bResponse );
    }
}
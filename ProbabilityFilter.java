/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.sf.jaer.eventprocessing.filter;

import net.sf.jaer.Description;
import net.sf.jaer.chip.AEChip;
import net.sf.jaer.event.EventPacket;
import net.sf.jaer.event.OutputEventIterator;
import net.sf.jaer.event.TypedEvent;
import net.sf.jaer.eventprocessing.EventFilter2D;
import net.sf.jaer.graphics.FrameAnnotater;

import com.jogamp.opengl.util.awt.TextRenderer;
import net.sf.jaer.DevelopmentStatus;
import net.sf.jaer.event.BasicEvent;
import net.sf.jaer.graphics.AEChipRenderer;
/**
 *
 * @author gary9555
 */
@Description("A probabilistic filter which renders the normal flow of moving edges")
public class ProbabilityFilter extends EventFilter2D{
    
    public ProbabilityFilter(AEChip chip) {
        super(chip);
        setPropertyTooltip("numDvsEventsToResetAccumulation", "sets number of dvs events to reset accumulation of image");
        setPropertyTooltip("showEventsAccumulatedBar", "shows a bar for num events accumulated");
        setPropertyTooltip("showTimeElapsedText", "shows text for time elapsed since last accumulation resst");
    }
    
    @Override
    public EventPacket<?> filterPacket(EventPacket<?> in) {
        
        return in;
    }
    
    @Override
    public void resetFilter() {
     
    }

    @Override
    public void initFilter() {
    }
    
    private float tmin = getFloat("tmin",10000f);
    private float tmax = getFloat("tmax",200000f);
    private float nres = getFloat("nres",20f);
    private float tdecay = getFloat("tdecay",10000f);
    private float tstd = getFloat("tstd",0.5f);
    
    ////////////////////////////////////////////////////////////////////
    ///////// ### Getters and setters ### //////////////////////////////
    ////////////////////////////////////////////////////////////////////    
    // getter for min time between spikes
    public float getTMin(){
      return tmin;
    }
    // getter for max time between spikes
    public float getTMax(){
      return tmax;
    }
    // getter for number of resolution
    public float getNRes(){
      return nres;
    }
    // getter for decay parameter
    public float getTDecay(){
      return tdecay;
    }
    // getter for standard deviation factor of the estimator
    public float getTStd(){
      return tstd;
    }

    /** setter for min time between spikes; updates GUI and saves as preference 
     @param NewFloat float value to set */
    public void setTMin(final float NewFloat) {
      putFloat("tmin",NewFloat);
      float OldValue = this.tmin;
      this.tmin = NewFloat;
      support.firePropertyChange("tmin",OldValue,NewFloat);
    }
    // setter for max time between spikes
    public void setTMax(final float NewFloat) {
      putFloat("tmax",NewFloat);
      float OldValue = this.tmax;
      this.tmax = NewFloat;
      support.firePropertyChange("tmax",OldValue,NewFloat);
    } 
    // setter for number of resolution 
    public void setNRes(final float NewFloat) {
      putFloat("nres",NewFloat);
      float OldValue = this.nres;
      this.nres = NewFloat;
      support.firePropertyChange("nres",OldValue,NewFloat);
    }
    // setter for decay parameter
    public void setTDecay(final float NewFloat) {
      putFloat("tdecay",NewFloat);
      float OldValue = this.tdecay;
      this.tdecay = NewFloat;
      support.firePropertyChange("tdecay",OldValue,NewFloat);
    }
    // setter for standard deviation factor of the estimator
    public void setTStd(final float NewFloat) {
      putFloat("tstd",NewFloat);
      float OldValue = this.tstd;
      this.tstd = NewFloat;
      support.firePropertyChange("tstd",OldValue,NewFloat);
    }
    ////////////////////////////////////////////////////////////////////
    ///////// ### End of getters and setters ### ///////////////////////
    ////////////////////////////////////////////////////////////////////
    
    
}


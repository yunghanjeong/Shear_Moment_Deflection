# Chassis Axle Deflection Calculation
## Shear-Moment Analysis and Deflection Visualizer

This repo is a refactored code of multilevel calculation performed during Chassis Axle Sensor prototyping during my time as Mechanical Engineer at PowerFleet, Inc (NASDAQ: PWFL).

![cooltractortrailor - getty images](https://static-25.sinclairstoryline.com/resources/media/c764120f-c6e0-48c6-ba00-0c8082020f10-large16x9_GettyImages822249792.jpg?1524838424831)

## Project Statement

It is often said that the trucks are the bloodline of America and as the world pushes for greener transportation many have pushed for more efficient methods of logistical transporation. To achieve higher efficiency set by many new legal guidelines the desire for more accurate tracking, measuring, and managing large tractor trailer has been fastly growing. This prototpye was envisioned to track the trailor's weight and load distribution and provide live statistics to the managing company via bluetooth and satelite communication. 

## Field Discoveries and Assumptions

* 40ft fixed bed, twin axle, dual wheel, is the most common type of trailer on the road.
* Each axle has their own leaf spring suspension. The suspension transfers the load from the frame to the axles. 
* There is a linear relationship between the trailor load and the suspension compression.
* There is measurable deflection of the axle when the trailer is under load

## Challenges

* The trailer is a statically indeterminate system. 
  * It can be simplified as a beam with a fixed support on one end (gooseneck or landing-gear) with 2 roller support on the other end (twin axle)
  * Shear-Moment analysis can be performed to calculate theoreitcal distribution of the load between all supports
* The trailer load is evenly distributed
  * Uneven distribution can be modeled as a combination of evenly distributed load and off-centered point load
* The load is evenly distributed between roadside and curbside (driver vs passenger side) suspension. 

## Methods

Solving engineering problems with numpy and matplotlib

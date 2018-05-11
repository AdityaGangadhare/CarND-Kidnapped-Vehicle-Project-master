/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

/*******************************************************************************
*                                 Initialization                               *
*******************************************************************************/
	// Initialize number of particles
	num_particles = 10;

	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	// This line creates a normal (Gaussian) distribution for y
	normal_distribution<double> dist_y(y, std_y);
	// This line creates a normal (Gaussian) distribution for theta
	normal_distribution<double> dist_theta(theta, std_theta);

/*******************************************************************************
*                                 Creating Particles                           *
*******************************************************************************/

	for (int i = 0; i < num_particles; i++) {

		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
		weights.push_back(1);
}

is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

/*******************************************************************************
*                          Creating normal distribution      		       *
*******************************************************************************/

	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(0, std_x);
	// This line creates a normal (Gaussian) distribution for y
	normal_distribution<double> dist_y(0, std_y);
	// This line creates a normal (Gaussian) distribution for theta
	normal_distribution<double> dist_theta(0, std_theta);

/*******************************************************************************
*                          Calculating state after a time step 		       *
*******************************************************************************/
	for (int i = 0; i < num_particles; i++) {

		if(fabs(yaw_rate) > 0.001){
			particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;
		}

		else {
			particles[i].x += velocity * delta_t * cos( particles[i].theta );
      		particles[i].y += velocity * delta_t * sin( particles[i].theta );
			particles[i].theta = particles[i].theta;
		}

		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
		particles[i].weight = 1;
		weights[i] = 1;
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {

/*******************************************************************************
*         Associating observed measurement to this particular landmark 	       *
*******************************************************************************/
	for (int i = 0; i < observations.size(); i++) { 

    double minDist;
    int id;

    for (int j = 0; j < predicted.size(); j++ ) { 

      double xDiff = observations[i].x - predicted[j].x;
      double yDiff = observations[i].y - predicted[j].y;

      double distance = xDiff * xDiff + yDiff * yDiff;

      // If the distance is less than min, store the id and update min.
      if ( distance < minDist || j==0 ) {
        minDist = distance;
        id = predicted[j].id;
      }
    }

    // Update the observation identifier.
    observations[i].id = id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	//   Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	//   NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i = 0; i < num_particles; i++) {

    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;

/*******************************************************************************
*  	                Find landmarks in particle's range	  	       *
*******************************************************************************/

    vector<LandmarkObs> landmarksWithinRange;
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      double dX = x - map_landmarks.landmark_list[j].x_f;
      double dY = y - map_landmarks.landmark_list[j].y_f;
      if ( dX*dX + dY*dY <= sensor_range*sensor_range ) {
        landmarksWithinRange.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f });
      }
    }


/*******************************************************************************
*     Convert observation coordinates from car coordinates to map coordinates  *
*******************************************************************************/
    vector<LandmarkObs> mObservations;
    for(unsigned int j = 0; j < observations.size(); j++) {
      double xx = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double yy = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      mObservations.push_back(LandmarkObs{ observations[j].id, xx, yy });
    }

    // Associate observation to landmark.
    dataAssociation(landmarksWithinRange, mObservations);


    // Reseting weight.
    particles[i].weight = 1.0;
    // Calculate new weights.
    for(unsigned int j = 0; j < mObservations.size(); j++) {
      double observationX = mObservations[j].x;
      double observationY = mObservations[j].y;
      int landmarkId = mObservations[j].id;
      double landmarkX, landmarkY;
      unsigned int k = 0;
      unsigned int nLandmarks = landmarksWithinRange.size();
      bool found = false;
      while( !found && k < nLandmarks ) {
        if ( landmarksWithinRange[k].id == landmarkId) {
          found = true;
          landmarkX = landmarksWithinRange[k].x;
          landmarkY = landmarksWithinRange[k].y;
        }
        k++;
      }

      double dX = observationX - landmarkX;
      double dY = observationY - landmarkY;
	//   std:: cout << observationX << endl;
	//   std:: cout << observationY << endl;
	//   std:: cout << landmarkX << endl;
	//   std:: cout << landmarkY << endl;
      double weight = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) * exp( -( dX*dX/(2*std_landmark[0]*std_landmark[0]) + (dY*dY/(2*std_landmark[1]*std_landmark[1])) ) );
      particles[i].weight *= weight;
	  weights[i] *= weight;
      
    }
  }
}

void ParticleFilter::resample() {
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

/*************************************************************************************
* Resample particles with replacement with probability proportional to their weight  *
**************************************************************************************/

	default_random_engine gen;
    discrete_distribution<int> distribution(weights.begin(),weights.end());    

    vector<Particle> resample_particles;

    for(int i=0;i<num_particles;i++)
    {
        resample_particles.push_back(particles[distribution(gen)]);
    }
    
    particles=resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

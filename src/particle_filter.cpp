/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <cassert>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
   if (is_initialized) 
   {
     return;
   }
  num_particles = 100;  // TODO: Set the number of particles
  
  default_random_engine gen;
  //creating normal distributions:
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  //
  for(int i = 0; i < num_particles; i++)
  {
    Particle new_particle;   
    new_particle.id = i;
    new_particle.x = dist_x(gen);
    new_particle.y = dist_y(gen);
    new_particle.theta = dist_theta(gen);
    new_particle.weight = 1.0;
    particles.push_back(new_particle);
    weights.push_back(new_particle.weight);    
  } 
  
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  default_random_engine gen;
    
  for(int i = 0; i < num_particles; i++)
  {
    double curr_theta = particles[i].theta;
    if(fabs(yaw_rate) < 0.0001)
    {
      particles[i].x += velocity * delta_t * cos(curr_theta);
      particles[i].y += velocity * delta_t * sin(curr_theta);
    }  
    else
    {
      particles[i].x += (velocity/yaw_rate)*(sin(curr_theta + yaw_rate*delta_t) - sin(curr_theta));
      particles[i].y += (velocity/yaw_rate)*(cos(curr_theta) - cos(curr_theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate*delta_t;    
    }
    
    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
  	normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
  	normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    // add random Gaussian noise:
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);;   
   
  }  
} 

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations,
                                    double sensor_range) 
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for(uint i=0; i<observations.size(); i++)
  {
    int landmark_id = -1;
    double min_dist = sensor_range*sensor_range;

    for(uint j=0; j<predicted.size(); j++)
    {
      double dist;
      dist = pow(observations[i].x - predicted[j].x,2) + pow(observations[i].y - predicted[j].y,2);
      if (dist < min_dist)
      {
        min_dist = dist;
        landmark_id = predicted[j].id;
      }    
    }
    observations[i].id = landmark_id;   
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
 /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */ 
  
  for(uint i=0; i<particles.size();i++)
  {
   vector<LandmarkObs> trans_observationss;
   for(uint j=0; j<observations.size(); j++)
   {
     // Step 1. Transformation to MAP coordiante system:
     LandmarkObs trans_observation;
     //trans_observation.id = j;
     //trans_observation.id = observations[j].id;
     trans_observation.id = -1;
     trans_observation.x = particles[i].x + cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y;
     trans_observation.y = particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
     trans_observationss.push_back(trans_observation);
   }
   //Step 2.Filter out landmarks that are out of sensor rahge of current particle
   vector<LandmarkObs> in_range_landmarks;
   for(uint k=0; k<map_landmarks.landmark_list.size(); k++)
   {
      double dist_land = pow(map_landmarks.landmark_list[k].x_f - particles[i].x,2);
      dist_land += pow(map_landmarks.landmark_list[k].y_f - particles[i].y,2);     
      if(dist_land < sensor_range*sensor_range)
      {
        LandmarkObs inrange_landmark;
        inrange_landmark.id = map_landmarks.landmark_list[k].id_i;
        inrange_landmark.x = map_landmarks.landmark_list[k].x_f;
        inrange_landmark.y = map_landmarks.landmark_list[k].y_f;
        in_range_landmarks.push_back(inrange_landmark);
      }     
   }
    
   // Step 3.Associate each measurement with a landmark identifier
   dataAssociation(in_range_landmarks, trans_observationss, sensor_range); // правильно ли работает
   // Step 4. Calculate final weight for particle
   //particles[i].weight = 1.0;
    
   double sigma_x = std_landmark[0];
   double sigma_y = std_landmark[1];
   double sigma_x2 = pow(sigma_x, 2);
   double sigma_y2 = pow(sigma_y, 2);
   double normalizer = 2*M_PI*sigma_x*sigma_y; 
   
   if(sigma_x2 == 0 || sigma_y2 == 0 || fabs(normalizer) < 0.0001)
   {
     normalizer = 0.0001;
   }
    
   double weight_meas = 1.0;
   
   for(uint n=0; n<trans_observationss.size(); n++)
    {         
       for(uint m=0; m<in_range_landmarks.size(); m++)
       {
         if(trans_observationss[n].id == in_range_landmarks[m].id) ///не выполняется это условие
         {
         	
           double x_dst2 = pow(trans_observationss[n].x - in_range_landmarks[m].x,2);
   		   double y_dst2 = pow(trans_observationss[n].y - in_range_landmarks[m].y,2);
           double exp_arg = x_dst2/(2.0*sigma_x2) + y_dst2/(2.0*sigma_y2);
         	
            weight_meas = weight_meas * (exp(-exp_arg))/normalizer;
         } // if        
       } // landmark m    
    } // observ n    
   
    particles[i].weight = weight_meas;
    weights[i] = weight_meas;
  } // particle i
}

void ParticleFilter::resample() 
{
  /*
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   vector<Particle> particles_res;
   default_random_engine gen;  
   uniform_int_distribution<int> particle_index(0, num_particles - 1);
   //generating index:
   int index = particle_index(gen);
   double beta = 0.0;
   double max_weight = numeric_limits<double>::min();
  
   for(uint i = 0; i < particles.size(); i++)
   {
     if(weights[i] > max_weight)
     {
       max_weight = weights[i];
     }
   }
  
   //uniform_real_distribution<double> rand_weight(0.0, max_weight);
  
  for(uint j=0; j<particles.size(); j++)
   {
     uniform_real_distribution<double> rand_weight(0.0, max_weight*2.0);
     beta += rand_weight(gen)*2.0;
     while(beta > weights[index])
     {
       beta -= weights[index];
       index = (index+1) % num_particles;      
     }  
     particles_res.push_back(particles[index]);
   }
  
   particles = particles_res; 
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
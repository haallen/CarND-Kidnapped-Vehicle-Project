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
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Create a normal (Gaussian) distribution for x, y, theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	num_particles = 50;

	for (int i = 0; i < num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Create a normal (Gaussian) distribution for x, y, theta
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	for (int i = 0; i < num_particles; i++) {
		if (fabs(yaw_rate)<0.001){
			particles[i].x += velocity*delta_t*cos(particles[i].theta);
		    particles[i].y += velocity*delta_t*sin(particles[i].theta);
		    particles[i].theta = particles[i].theta;
		}else{
			particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
			particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;
		}

		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i<observations.size(); i++){
		double min_dist = numeric_limits<double>::max();
		int closest_idx = -1;
		for (int j = 0; j<predicted.size(); j++){
			double calc_dist = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
			if (calc_dist < min_dist){
				min_dist = calc_dist;
				closest_idx = j;//predicted[j].id;
			}
		}
		observations[i].id = closest_idx;
		//std::cout<<"obs "<<observations[i].x<<" "<<observations[i].y<< " "<<observations[i].id<<std::endl;
		//std::cout<<"pred "<<predicted[closest_idx].x<<" "<<predicted[closest_idx].y<< " "<<predicted[closest_idx].id<<std::endl;
		//std::cout<<" closest idx "<< closest_idx<<std::endl;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	for (int j = 0; j<num_particles; j++){

		//find the landmarks that are within sensor range
		std::vector<LandmarkObs> range_landmarks;
		for (int i = 0; i<map_landmarks.landmark_list.size();i++){
			if(dist(map_landmarks.landmark_list[i].x_f,map_landmarks.landmark_list[i].y_f,particles[j].x,particles[j].y)<sensor_range){
				range_landmarks.push_back(LandmarkObs{map_landmarks.landmark_list[i].id_i,map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f});
			}
		}

		//rotate the observations to map coordinate frame
		std::vector<LandmarkObs> rot_observations;
		for (int i = 0; i<observations.size();i++){
			double rot_x = particles[i].x + cos(particles[j].theta)*observations[i].x - sin(particles[j].theta)*observations[i].y;
			double rot_y = particles[i].y + sin(particles[j].theta)*observations[i].x + cos(particles[j].theta)*observations[i].y;
			int rot_id = observations[i].id;
			rot_observations.push_back(LandmarkObs{rot_id, rot_x, rot_y});
		}

		//associate each observation to a landmark
		dataAssociation(range_landmarks,rot_observations);

		//calculate normalization term
		double gauss_norm= (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));

		particles[j].weight = 1.0;

		for (int i = 0; i<rot_observations.size();i++){

			// calculate exponent
			double mu_x = range_landmarks[rot_observations[i].id].x;
			double mu_y = range_landmarks[rot_observations[i].id].y;
			double diff_x = rot_observations[i].x - mu_x;
			double diff_y = rot_observations[i].y - mu_y;
			double exponent1 = (diff_x*diff_x)/(2*std_landmark[0]*std_landmark[0]);
			double exponent2 = (diff_y*diff_y)/(2*std_landmark[1]*std_landmark[1]);
			double exponent = exponent1+exponent2;
			//calculate weight using normalization terms and exponent
			particles[j].weight *= (gauss_norm *exp(-exponent));
		}
		//std::cout<<"particle "<<j<<" "<<particles[j].weight<<std::endl;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	uniform_int_distribution<> distribution_i(0,num_particles-1);
	int index = distribution_i(gen);

	vector<Particle> sampled_particles;
	double beta = 0.0;

	double mw = 0.0;
	for (int j = 0; j<num_particles; ++j){
		if (particles[j].weight > mw){
			mw = particles[j].weight;
		}
	}
	uniform_real_distribution<int> distribution_r(0.0,mw);

	for (int j = 0; j<num_particles; j++){
		beta += distribution_r(gen)*2.0;
		while (beta > particles[index].weight){
			beta -= particles[index].weight;
			index += 1;
		}
		sampled_particles.push_back(particles[index]);
	}
	particles = sampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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

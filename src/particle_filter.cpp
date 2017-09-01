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
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	default_random_engine gen;

	for (int i = 0; i < num_particles; i++) {

		Particle sample;
		sample.id=i;
		sample.x=dist_x(gen);
		sample.y=dist_y(gen);
		sample.theta=dist_theta(gen);
		sample.weight=1.0;

		particles.push_back(sample);

	}

	is_initialized=true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // define normal distributions for sensor noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
	default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {

    if (fabs(yaw_rate) < 0.000001) {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }
    else {
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate * delta_t;
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

	int min_index=0;
	double distance= 0;
	for (int i = 0; i < observations.size(); i++) {
		double min=999999;
		for(int j=0; j< predicted.size();j++){
				distance=dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
				if(distance<min){
					min=distance;
					min_index=j;
				}
		}
		observations[i].id=predicted[min_index].id;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html



	for (int i = 0; i < num_particles; i++) {

		std::vector<LandmarkObs> predicted;
		std::vector<LandmarkObs> t_observations;

		for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			if(dist(particles[i].x,particles[i].y,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f)<sensor_range){
				LandmarkObs predicted_landmark;
				predicted_landmark.id=map_landmarks.landmark_list[j].id_i;
				predicted_landmark.x=map_landmarks.landmark_list[j].x_f;
				predicted_landmark.y=map_landmarks.landmark_list[j].y_f;
				predicted.push_back(predicted_landmark);

			}
		}
		for (int j = 0; j < observations.size(); j++) {
			LandmarkObs transformed_ob;
			transformed_ob.id=observations[j].id;
			transformed_ob.x= particles[i].x+ cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
			transformed_ob.y= particles[i].y+ sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
			t_observations.push_back(transformed_ob);
		}

		dataAssociation(predicted, t_observations);

		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		for (int j = 0; j < t_observations.size(); j++) {
			associations.push_back(t_observations[j].id);
			sense_x.push_back(t_observations[j].x);
			sense_x.push_back(t_observations[j].y);
		}
		SetAssociations(particles[i], associations, sense_x, sense_y);
		particles[i].weight=1.0;
		//cout<<particles[i].weight<<endl;
		for(int j=0;j<t_observations.size();j++){
			int assoc_index=0;
			for(int k=0;k<predicted.size();k++){
				if(predicted[k].id==t_observations[j].id){
					assoc_index=k;
					break;
				}
			}

			particles[i].weight = particles[i].weight *( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) * exp( -( pow(predicted[assoc_index].x-t_observations[j].x,2)/(2*pow(std_landmark[0], 2)) + (pow(predicted[assoc_index].y-t_observations[j].y,2)/(2*pow(std_landmark[1], 2))) ) );
			//cout<<"prx"<<predicted[assoc_index].x<<endl;
			cout<<endl;
			cout<<"pid"<<predicted[assoc_index].id<<endl;
			cout<<"prx"<<predicted[assoc_index].x<<endl;
			cout<<"pry"<<predicted[assoc_index].y<<endl;
			cout<<"ox"<<t_observations[j].x<<endl;
			cout<<"ox"<<t_observations[j].x<<endl;
			cout<<"oy"<<t_observations[j].y<<endl;
		}
		cout<<endl;
		//cout<<particles[i].weight<<endl;
		weights.push_back(particles[i].weight);

	}
}
void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::vector<Particle> new_particles;
	default_random_engine gen;
	std::discrete_distribution<> distribution(weights.begin(),weights.end());

	for(int i=0;i<num_particles;i++){
		new_particles.push_back(particles[distribution(gen)]);
	}
	particles=new_particles;
	weights.clear();
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

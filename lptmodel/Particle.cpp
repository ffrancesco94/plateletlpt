#include "Particle.h"
#include "io.h"
#include "macros.h"
#include "Fluid.h"
#include "typedefs.h"
#include "DynamicFactory.hh"

// Explicit instantiation of factory
//template class DynamicFactory<Particle>;

ParticleFactory & particleFactory()
{
	static ParticleFactory particleFactory;
	return particleFactory;
}

// Particle
Particle::Particle(const Particle & rhs)
// : density_(rhs.density()), radius_(rhs.radius())
{
	// Clone particle forces
	density_ = rhs.density();
	radius_ = rhs.radius();
	for(auto && particleForce : rhs.particleForces())
		this->particleForces().emplace_back(std::unique_ptr<ParticleForce>(particleForce->clone()));
}

Particle & Particle::operator=(const Particle & rhs)
{
	if(this != &rhs) {
		this->particleForces().clear();
		for(auto && particleForce : rhs.particleForces())
			this->particleForces().emplace_back(std::unique_ptr<ParticleForce>(particleForce->clone()));
		radius() = rhs.radius();
		density() = rhs.density();
	}
	return *this;
}
void Particle::writeBinary(std::ostream & os) const
{
	write_to_stream(os, typeId());
	write_to_stream(os, id_);
	write_to_stream(os, position());
	write_to_stream(os, velocity());
	write_to_stream(os, shear_);
	write_to_stream(os, pas_);
	write_to_stream(os, dose_);
	write_to_stream(os, age_);
	write_to_stream(os, isAlive_);
	write_to_stream(os, collisionCount_);
	write_to_stream(os, injectionTime_);
}

void Particle::readBinary(std::istream & is)
{
	read_from_stream(is, id_);
	read_from_stream(is, position());
	read_from_stream(is, velocity());
	read_from_stream(is, shear_);
	read_from_stream(is, pas_);
	read_from_stream(is, dose_);
	read_from_stream(is, age_);
	read_from_stream(is, isAlive_);
	read_from_stream(is, collisionCount_);
	read_from_stream(is, injectionTime_);
}

// Tracer particle
int TracerParticle::typeId_ = registerParticleTypeToFactory<TracerParticle>("TracerParticle");
std::string TracerParticle::typeName_ = "TracerParticle";

TracerParticle * TracerParticle::clone() const
{ 
	return new TracerParticle(*this); 
}

void TracerParticle::updateMomentum(scalar dt, const Vector & fluidVelocity, const Matrix & shear, const Fluid & fluid, const Vector & fluidVorticity)
{
	this->velocity() = fluidVelocity;
}


// Material particle
int MaterialParticle::typeId_ = registerParticleTypeToFactory<MaterialParticle>("MaterialParticle");
std::string MaterialParticle::typeName_ = "MaterialParticle";

// MaterialParticle::MaterialParticle(const MaterialParticle & rhs)
// // : density_(rhs.density()), radius_(rhs.radius())
// {
// 	// Clone particle forces
// 	density_ = rhs.density();
// 	radius_ = rhs.radius();
// 	for(auto && particleForce : rhs.particleForces())
// 		this->particleForces().emplace_back(std::unique_ptr<ParticleForce>(particleForce->clone()));
// }

// MaterialParticle & MaterialParticle::operator=(const MaterialParticle & rhs)
// {
// 	if(this != &rhs) {
// 		this->particleForces().clear();
// 		for(auto && particleForce : rhs.particleForces())
// 			this->particleForces().emplace_back(std::unique_ptr<ParticleForce>(particleForce->clone()));
// 		radius() = rhs.radius();
// 		density() = rhs.density();
// 	}
// 	return *this;
// }

MaterialParticle * MaterialParticle::clone() const
{ 
	return new MaterialParticle(*this); 
}

void MaterialParticle::updateMomentum(scalar dt, const Vector & fluidVelocity, const Matrix & shear, const Fluid & fluid, const Vector & fluidVorticity)
{
	// Compute total force
	Vector totalForce(0, 0, 0);
	// New stuff
	const int drag_id = StokesDrag().typeId();
	const int saffman_id = SaffmanForce().typeId();
	const scalar mass = 4. / 3. * M_PI * radius() * radius() * radius() * density();
	Vector alpha(0,0,0);
	scalar beta = 0.;
	// End new
	ParticleForceData forceData{fluidVelocity, velocity(), shear, fluid, radius(), density(), position(), fluidVorticity};
	// std::cout << "Fluid velocity: " << fluidVelocity << std::endl;
	// std::cout << "Particle velocity: " << velocity() << std::endl;
	// std::cout << "Radius: " << radius() << std::endl;
	// std::cout << "Computing " << this->particleForces().size() << " forces" << std::endl;
	for(auto && particleForce : this->particleForces()) {
		totalForce += particleForce->getParticleForce(forceData);
		// New stuff
		if (particleForce->typeId() == drag_id) {
			// std::cout << "Drag : " << particleForce->getParticleForce(forceData).norm() << std::endl;
			// std::cout << "Slip velocity: " << (forceData.fluidVelocity - velocity()).norm() << std::endl;
			if ((forceData.fluidVelocity - velocity()).norm() > 1e-15) {
				beta = particleForce->getParticleForce(forceData).norm() / (forceData.fluidVelocity - velocity()).norm();
			}
			alpha += beta * forceData.fluidVelocity;
		} else if (particleForce->typeId() == saffman_id) {
			alpha += particleForce->getParticleForce(forceData);
		}
		// End new
	}
	// New stuff
	alpha /= mass;
	beta /= mass;
	// I want to write the equation of motion as du/dt = (\alpha - \beta*u).
	// That way I can write du = (\alpha - \beta*u)dt/(1+\beta*dt). See what OF does!
	// In this case, \alpha = (D/u_s*u_f + L)/m_p
	//				 \beta = D/(m_p*u_s)
	const Vector du = (alpha - beta * this->velocity()) * dt / (1. + beta * dt);
	this->velocity() += du;
	// End new
	// Explicit Euler velocity
	//this->velocity() += 3. * dt * totalForce / (4. * density() * radius() * radius() * radius() * M_PI);

	// this->velocity() = fluidVelocity;
}

int MaterialParticle::typeId() const
{ 
	return typeId_; 
}

void MaterialParticle::fromJSON(const json & jsonObject)
{
	Particle::fromJSON(jsonObject);
	density() = jsonObject.at("density");
	radius() = jsonObject.at("radius");
	particleForces() = particleForceFactory().createVectorFromJSON(jsonObject.at("forces"));
}

void MaterialParticle::writeBinary(std::ostream & out) const
{
	Particle::writeBinary(out);
	write_to_stream(out, density());
	write_to_stream(out, radius());

	// Write particle forces
	write_to_stream<int>(out, particleForces().size());
	for(auto && particleForce : particleForces())
		particleForce->writeBinary(out);
}

void MaterialParticle::readBinary(std::istream & in)
{
	Particle::readBinary(in);
	read_from_stream(in, density());
	read_from_stream(in, radius());

	int numberOfParticleForces;
	particleForces().clear();
	read_from_stream(in, numberOfParticleForces);
	for(int i = 0; i < numberOfParticleForces; ++i)
		particleForces().emplace_back(particleForceFactory().createFromStream(in));
}

// No momentum update particle
int NoMomentumUpdateParticle::typeId_ = registerParticleTypeToFactory<NoMomentumUpdateParticle>("NoMomentumUpdateParticle");
std::string NoMomentumUpdateParticle::typeName_ = "NoMomentumUpdateParticle";
	

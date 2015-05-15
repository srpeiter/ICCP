//============================================================================
// Name        : renderingtest.cpp
// Author      : Sarwan Peiter
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Rendering.h"
//float *pos[3];
//const int N=40000;

GLuint m_texStar;
GLFWwindow* window;
//int m_idxSnapshot;

static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

/*

void SaveToTGA(int idx)
{
  if (idx==-1)
    m_idxSnapshot++;
  else
    m_idxSnapshot = idx;

  std::stringstream ss;
  ss << "frame_" << std::setw(5) << std::setfill('0') << m_idxSnapshot << ".tga";
  SaveToTGA2(ss.str());
}

void SaveToTGA2(const std::string &sName)
{
	int width=0, heigth=0;
  using std::ios;

  glViewport(0, 0, width, heigth);

  int nSize = width * heigth * 3;

  GLubyte pixels[nSize];
  glReadPixels(0, 0, width, heigth, GL_RGB, GL_UNSIGNED_BYTE, pixels);

  std::string sFile;
  if (sName.length())
    sFile = sName;
  else
  {
    // use default name with time stamp
    time_t t = time(NULL);
    struct tm *tmp = localtime(&t);
    if (tmp==NULL)
      sFile = "snapshot.tga";
    else
    {
      char szTime[1024];
      if (strftime(szTime, sizeof(szTime), "snapshot_%Y%m%d_%H%M%S.tga", tmp) == 0)
        sFile = "snapshot.tga";
      else
        sFile = szTime;
    }
  }

  std::fstream file(sFile.c_str(), ios::out|ios::binary|ios::trunc);
  char TGAheader[12] = { 0,0,2,0,0,0,0,0,0,0,0,0 };
  char header[6] = { width  % 256,
                     width  / 256,
                     heigth % 256,
                     heigth / 256,
                     24,
                     0 };
  file.write(TGAheader, sizeof(TGAheader));
  file.write(header, sizeof(header));

  //convert to BGR format
  for (int i=0; i<nSize; i+=3)
    std::swap(pixels[i], pixels[i+2]);

  file.write(reinterpret_cast<char*>(pixels), nSize);
  file.close();
}

*/

void initGLFW()
{
	glfwSetErrorCallback(error_callback);
	    if (!glfwInit())
	        exit(EXIT_FAILURE);
	   // window = glfwCreateWindow(2560, 1440, "Colliding Galaxy Simulation", NULL, NULL);
	    window = glfwCreateWindow(640, 480, "Colliding Galaxy Simulation", NULL, NULL);
	    if (!window)
	    {
	        glfwTerminate();
	        exit(EXIT_FAILURE);
	    }

	    // make current window OPENGL context or else OPENGL will not work
	    glfwMakeContextCurrent(window);

	    // time interval at which the render buffer(back buffer) is brought to the front buffer
	    glfwSwapInterval(1);

	    //terminate when escape-key is pressed
	    glfwSetKeyCallback(window, key_callback);
}

/* Initialize OpenGL Graphics */
void initGL() {


	glEnable (GL_LINE_SMOOTH);
	  glEnable (GL_BLEND);
	  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	  glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	  glLineWidth (1);

	  glClear(GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT);
	  glEnable(GL_POINT_SPRITE);
	  glDisable(GL_DEPTH_TEST);
	  glClearColor(0.0f, 0.0f, 0.03f, 0.0f);  // black background
	  //SetCameraOrientation(Vec3D(0, 1, 0));

	  glMatrixMode(GL_MODELVIEW);
	  glLoadIdentity();

}

void initPointSpriteExt()
{

	SDL_Surface *tex;

	   // texture loading taken from
	// http://gpwiki.org/index.php/SDL:Tutorials:Using_SDL_with_OpenGL
	 // tex = SDL_LoadBMP("particle.bmp");
	 tex = SDL_LoadBMP("Star.bmp");
	  if (!tex)
	    throw std::runtime_error("can't load star texture (particle.bmp).");

		// Check that the image's width is a power of 2
		if (tex->w & (tex->w - 1))
			throw std::runtime_error("texture width is not a power of 2.");

		// Also check if the height is a power of 2
		if (tex->h & (tex->h - 1))
			throw std::runtime_error("texture height is not a power of 2.");

	  // get the number of channels in the SDL surface
	  GLint  nOfColors = tex->format->BytesPerPixel;
	  GLenum texture_format;
	  if (nOfColors == 4)     // contains an alpha channel
	  {
	    if ( tex->format->Rmask == 0x000000ff)
	      texture_format = GL_RGBA;
	    else
	      texture_format = GL_BGRA;
	  }
	  else if (nOfColors == 3)     // no alpha channel
	  {
	    if ( tex->format->Rmask == 0x000000ff)
	      texture_format = GL_RGB;
	    else
	      texture_format = GL_BGR;
	   }
	   else
	     throw std::runtime_error("image is not truecolor");

	  // Have OpenGL generate a texture object handle for us
		glGenTextures(1, &m_texStar);

	  // Bind the texture object
		glBindTexture( GL_TEXTURE_2D, m_texStar );

		// Set the texture's stretching properties
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		// Edit the texture object's image data using the information SDL_Surface gives us
		glTexImage2D(GL_TEXTURE_2D,
		             0,
		             nOfColors,
		             tex->w,
		             tex->h,
		             0,
	               texture_format,
	               GL_UNSIGNED_BYTE,
	               tex->pixels );

}

void Pre_Render()
{
	initGLFW();
	initGL();
	initPointSpriteExt();
}

void Render(std::vector<prop> coord, int num)
{

	glClear(GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT);

	int width, height;
	float r = 1;
	float g = 1;
	float b = 1;



	  // start rendering

	    float maxSize = 0.0f;

	    //------------------------------rendering------------------






	        glfwGetFramebufferSize(window, &width, &height); // not important

	        glViewport(0, 0, width, height);  // view settings




            glMatrixMode( GL_MODELVIEW );

	        glLoadIdentity();

	        glScalef(0.04, 0.04, 0.04);


	        glViewport(0, 0, width, height);

	   //     std::cout << width << " " ;
	     //   glMatrixMode(GL_PROJECTION);
	       // glLoadIdentity();
	        //glOrtho(-width/height, width/height, -1.f, 1.f, 1.f, -1.f);
	        //glOrtho(-1, 0, -0.f, 0.f, 0.f, -0.f);
	        //glMatrixMode(GL_MODELVIEW);
	        //glLoadIdentity();
	        glRotatef((float) glfwGetTime() * 5.f, 2.f, 4.f, 0.6f);



	        glBindTexture(GL_TEXTURE_2D, m_texStar);

	         glGetFloatv( GL_POINT_SIZE_MAX_ARB, &maxSize );

	         glTexEnvf(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);

	         glEnable(GL_POINT_SPRITE_ARB);
	         glEnable(GL_TEXTURE_2D);       // point sprite texture support
	         glEnable(GL_BLEND);            // soft blending of point sprites
	         glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	          // Render a color-cube consisting of 6 quads with different colors





	       // glMatrixMode(GL_MODELVIEW);


	          glGetFloatv( GL_POINT_SIZE_MAX_ARB, &maxSize );

	          glTexEnvf(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);

	          glEnable(GL_POINT_SPRITE_ARB);
	          glEnable(GL_TEXTURE_2D);       // point sprite texture support
	          glEnable(GL_BLEND);            // soft blending of point sprites
	          glBlendFunc(GL_SRC_ALPHA, GL_ONE);



	          glPointSize(4);
	            glBegin(GL_POINTS);



	          for (int i=0; i<num; ++i)
	          {
	          glPushMatrix();

	          glTranslatef(coord[i].x,coord[i].y, coord[i].z);



	           //glColor4f(1,coord[i].mass*g*10000,coord[i].mass*b*10000, 0.6);
	         glColor4f(coord[i].xcol,coord[i].ycol,coord[i].zcol, 06);

	            glVertex3f(coord[i].x,coord[i].y, coord[i].z);

	            glPopMatrix();
	          }
	          glEnd();


	          glDisable(GL_POINT_SPRITE_ARB);
	          glDisable(GL_BLEND);
	          glDisable(GL_TEXTURE_2D);

	        glfwSwapBuffers(window);
	        glfwPollEvents();
	    //}


}

void Post_Render()
{
	glDeleteTextures(1,&m_texStar);
	glfwDestroyWindow(window);
    glfwTerminate();
}

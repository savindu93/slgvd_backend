from django.contrib import admin
from django.urls import path, include
from api.views import CustomTokenObtainPairView
from rest_framework import permissions
from rest_framework_simplejwt.views import TokenRefreshView
from drf_yasg.views import get_schema_view
from drf_spectacular.views import SpectacularAPIView, SpectacularRedocView, SpectacularSwaggerView
from drf_yasg import openapi


urlpatterns = [
    path('admin/', admin.site.urls),
    path('api-auth/', include('rest_framework.urls')),
    path('api/token/', CustomTokenObtainPairView.as_view()),
    path('api/token/refresh/', TokenRefreshView.as_view()),
    path('api/', include('api.urls')),
    path('schema/', SpectacularAPIView.as_view(), name = 'schema'),
    path('swagger/', SpectacularSwaggerView.as_view(url_name = 'schema')),
    path('redoc/', SpectacularRedocView.as_view(url_name = 'schema')),
]

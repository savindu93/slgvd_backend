from django.urls import path
from api.consumer import ProgressConsumer

websocket_urlpatterns = [
    path('ws/progress/', ProgressConsumer.as_asgi()),
]
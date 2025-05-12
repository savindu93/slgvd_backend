from django.contrib.auth.backends import ModelBackend
from django.contrib.auth import get_user_model

# EmailBackend() redefines the authenticate() during login where instead of using the username as a field for authentication it uses the email of the user.
class EmailBackend(ModelBackend):

    def authenticate(self, request, email, password):

        User = get_user_model()

        try:
            user = User.objects.get(email=email)
                    
        except User.DoesNotExist:
            return None
        
        if user.password == password and self.user_can_authenticate(user):
            return user

        
